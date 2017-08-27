package Package;
//
import java.util.ArrayList;
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.basic.InterpolationModel;
import miniufo.basic.InterpolationModel.Type;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Typhoon;

//
public class LifeCycle{
	//
	private static int ystr=1989;
	private static int yend=2009;
	
	//private static int effectiveCount=0;
	private static int effectiveTC   =0;
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String tranges="time=1Jan"+ystr+"-31Dec"+yend;
	
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		List<float[]> winds=new ArrayList<float[]>();
		List<float[]> press=new ArrayList<float[]>();
		List<float[]> ssts =new ArrayList<float[]>();
		List<float[]> vwss =new ArrayList<float[]>();
		List<float[]> efcss=new ArrayList<float[]>();
		List<float[]> efcls=new ArrayList<float[]>();
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,25,2,0.3f,-650));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v","sst");
			Variable[] shrs=dm.cVerticalWindShear(vars[0],vars[1]);
			
			Variable shrsum=dm.cRadialAverage(shrs[0],1,15).anomalizeX();
			Variable shrsvm=dm.cRadialAverage(shrs[1],1,15).anomalizeX();
			
			Variable vwsm=shrsum.hypotenuse(shrsvm);
			Variable sstm=dm.cRadialAverage(vars[2],1,15).anomalizeX().minusEq(273.15f);
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efclm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),15,24);	// 500-800 km
			Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]), 9,18);	// 300-600 km
			
			boolean[] wind=IntensityModel.greaterEqualThan(tr.getWinds(),17.2f);
			
			if(IntensityModel.getCycles(wind)!=1)
				System.out.println(String.format("%7s (%4s)",tr.getName(),tr.getID()));
			else{
				int[] tag=extractValidDim(tr.getWinds(),17.2f);
				
				float[] tmp=extractWinds(tr.getWinds(),tag);
				if(tmp.length>1){
					effectiveTC++;
					winds.add(tmp);
					press.add(extractPres(tr.getPressures(),tag));
					 ssts.add(extractPres(sstm.getData()[0][0][0],tag));
					 vwss.add(extractPres(vwsm.getData()[0][0][0],tag));
					efcss.add(extractPres(efcsm.getData()[1][0][0],tag));
					efcls.add(extractPres(efclm.getData()[1][0][0],tag));
				}
			}
		}
		
		float[] comWind=stdComposite(winds);
		float[] comPres=stdComposite(press);
		float[] comSST =stdComposite(ssts );
		float[] comVWS =stdComposite(vwss );
		float[] comEFCS=stdComposite(efcss);
		float[] comEFCL=stdComposite(efcls);
		
		System.out.println("\nComposite length: "+comWind.length+"\n");
		
		for(int l=0;l<comWind.length;l++)
		System.out.println(String.format(
			"%6.2f  %7.2f  %4.1f  %4.1f  %5.1f  %5.1f",
			comWind[l],comPres[l],comSST[l],comVWS[l],comEFCS[l],comEFCL[l]
		));
		
		System.out.println();
		System.out.println("effective TCs  :"+effectiveTC);
	}
	
	private static float[] stdComposite(List<float[]> ls){
		int lenMax=0;
		for(float[] f:ls) if(f.length>lenMax) lenMax=f.length;
		
		float[] re=new float[lenMax];
		
		for(float[] f:ls)
		if(f.length<lenMax){
			float[] tmp=InterpolationModel.interp1D(f,lenMax,Type.LINEAR);
			
			for(int l=0;l<lenMax;l++) re[l]+=tmp[l];
			
		}else{
			for(int l=0;l<lenMax;l++) re[l]+=f[l];
		}
		
		if(ls.size()!=effectiveTC) throw new IllegalArgumentException();
		
		for(int l=0;l<lenMax;l++) re[l]/=ls.size();
		
		return re;
	}
	
	private static float[] extractWinds(float[] wind,int[] tag){
		float[] re=new float[tag[1]];
		
		for(int l=0;l<tag[1];l++) re[l]=wind[tag[0]+l];
		
		return re;
	}
	
	private static float[] extractPres(float[] pres,int[] tag){
		float[] re=new float[tag[1]];
		
		for(int l=0;l<tag[1];l++) re[l]=pres[tag[0]+l];
		
		return re;
	}
	
	// [0] is start tag, [1] is length
	private static int[] extractValidDim(float[] wind,float threshold){
		int[] re=new int[2];
		
		for(int l=0;l<wind.length;l++)
		if(wind[l]>=threshold){ re[0]=l;break;}
		
		for(int l=0;l<wind.length;l++)
		if(wind[l]>=threshold) re[1]++;
		
		return re;
	}
}
