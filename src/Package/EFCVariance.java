package Package;
//
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.basic.ArrayUtil;
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
public class EFCVariance{
	//
	private static int effectiveCount=0;
	private static int effectiveTC   =0;
	
	private static float minLat= 90;
	private static float minLon=360;
	private static float maxLat=-90;
	private static float maxLon=  0;
	
	private static final boolean noDepress=false;	// wind > 17.2 m/s
	private static final boolean noLanding=false;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	
	private static final DataSets dsets=DataSets.JMA;
	//private static final String respath="d:/Data/PhD/Climatology/"+dsets+"/";
	private static final String tranges="time=1Jan1989-31Dec2009";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		
		double sumEFCS=0;
		double sumFFCS=0;
		double sumEFCL=0;
		double sumFFCL=0;
		
		int count=0;
		int cc=0;
		for(Typhoon tr:ls){
			if(++cc%20==0) System.out.print(".");
			
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,28,2,0.3f,-650));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			float tmp;
			tmp=ArrayUtil.getMin(dd.getLat());	if(tmp<minLat) minLat=tmp;
			tmp=ArrayUtil.getMax(dd.getLat());	if(tmp>maxLat) maxLat=tmp;
			tmp=ArrayUtil.getMin(dd.getLon());	if(tmp<minLon) minLon=tmp;
			tmp=ArrayUtil.getMax(dd.getLon());	if(tmp>maxLon) maxLon=tmp;
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v");
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getZonalVelocity(),tr.getMeridionalVelocity(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]), 9,18);	// 300-600 km
			Variable efclm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),18,27);	// 600-900 km
			
			//ffcsm.multiplyEq(1e9f);
			//ffclm.multiplyEq(1e9f);
			
			float[] EFCS=efcsm.getData()[1][0][0];
			float[] EFCL=efclm.getData()[1][0][0];
			
			double tmpEFCS=0;
			double tmpFFCS=0;
			double tmpEFCL=0;
			double tmpFFCL=0;
			for(int i=0;i<EFCS.length;i++){
				double
				v=EFCS[i];	tmpEFCS+=v*v;
				v=EFCL[i];	tmpEFCL+=v*v;
			}
			
			sumEFCS+=tmpEFCS;
			sumFFCS+=tmpFFCS;
			sumEFCL+=tmpEFCL;
			sumFFCL+=tmpFFCL;
			count+=EFCS.length;
			effectiveTC++;
		}
		
		printResult(effectiveTC,effectiveCount);
		
		System.out.println("sum of EFCS: "+sumEFCS+" ("+Math.sqrt(sumEFCS/(count-1))+")");
		System.out.println("sum of FFCS: "+sumFFCS+" ("+Math.sqrt(sumFFCS/(count-1))+")");
		System.out.println("sum of EFCL: "+sumEFCL+" ("+Math.sqrt(sumEFCL/(count-1))+")");
		System.out.println("sum of FFCL: "+sumFFCL+" ("+Math.sqrt(sumFFCL/(count-1))+")");
		System.out.println("total count: "+count);
		
		//System.out.println(String.format(
		//	"\nsum of EFC: %.4f (%.4f)\nsum of FFC: %.4f (%.4f)\ntotal count: %d",
		//	sumEFC,sumEFC/count,sumFFC,sumFFC/count,count
		//));
	}
	
	static void printResult(int effectiveTC,int effectiveCount){
		System.out.println("\n\nwithin the region:\nlat["+
			String.format("%6.2f",minLat)+" N, "+
			String.format("%6.2f",maxLat)+" N]\nlon["+
			String.format("%6.2f",minLon)+" E, "+
			String.format("%6.2f",maxLon)+" E]\n"
		);
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
	}
}
