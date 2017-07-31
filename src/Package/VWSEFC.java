package Package;
//
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.lagrangian.Typhoon;


//
public class VWSEFC{
	//
	private static int VWSLPEFCSCount=0;
	private static int VWSSPEFCSCount=0;
	
	private static final float EFCThreshold=20f;
	private static final float VWSThreshold=10f;
	
	private static final boolean noDepress=true;	// wind > 17.2 m/s
	private static final boolean noLanding=true;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Climatology/"+dsets+"/";
	private static final String tranges="time=1Jan1989-31Dec2009";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		Script VWSLPEFCS=new Script(ctl,respath+"VWSLPEFCS.gs");
		Script VWSSPEFCS=new Script(ctl,respath+"VWSSPEFCS.gs");
		
		int cc=0;
		for(Typhoon tr:ls){
			if(++cc%20==0) System.out.print(".");
			
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
			dm.cStormRelativeAziRadVelocity(tr.getZonalVelocity(),tr.getMeridionalVelocity(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]), 9,18);	// 300-600 km
			
			ct=new CoordinateTransformation(new SphericalSpatialModel(ctl2),csm);
			Variable lsmm=dm.cRadialAverage(ct.transToCylindricalInvariantly(lsm),1,6).anomalizeX();
			
			boolean[] wind=noDepress?	// wind >= 17.2 m/s
				IntensityModel.greaterEqualThan(tr.getWinds(),17.2f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] lsmb=noLanding?	// no land within 200 km
				IntensityModel.lessThan(lsmm.getData()[0][0][0],1e-9f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] llon=withinWNP?	// lons >= 100E
				IntensityModel.greaterEqualThan(tr.getLongitudes(),100):IntensityModel.newBooleans(tr.getTCount());
			boolean[] rlon=withinWNP?	// lons <= 190E
				IntensityModel.lessEqualThan(tr.getLongitudes(),190):IntensityModel.newBooleans(tr.getTCount());
			
			boolean[] PEFCS=IntensityModel.greaterEqualThan(efcsm.getData()[1][0][0],EFCThreshold);
			boolean[] VWSL =IntensityModel.greaterEqualThan(vwsm.getData()[0][0][0],VWSThreshold);
			boolean[] VWSS =IntensityModel.lessThan(vwsm.getData()[0][0][0],VWSThreshold);
			
			boolean[] combVWSLPEFCS=IntensityModel.combination(wind,VWSL,PEFCS,lsmb,llon,rlon);
			boolean[] combVWSSPEFCS=IntensityModel.combination(wind,VWSS,PEFCS,lsmb,llon,rlon);
			
			float[] EFCS=efcsm.getData()[1][0][0];
			float[] VWSM= vwsm.getData()[0][0][0];
			float[] SSTM= sstm.getData()[0][0][0];
			float[] pote=IntensityModel.cPotential(sstm.getData()[0][0][0],tr.getStormRelativeWinds());
			
			for(int l=0;l<tr.getTCount();l++){
				if(combVWSLPEFCS[l]){ VWSLPEFCS.addDrawScript(tr,l,EFCS,VWSM,SSTM,pote);VWSLPEFCSCount++;}
				if(combVWSSPEFCS[l]){ VWSSPEFCS.addDrawScript(tr,l,EFCS,VWSM,SSTM,pote);VWSSPEFCSCount++;}
			}
		}
		
		VWSLPEFCS.flush();	VWSSPEFCS.flush();
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("\nthreshold: "+EFCThreshold+" m s^-1 / day");
		System.out.println("effective PEFCS count:"+VWSLPEFCSCount);
		System.out.println("effective NEFCS count:"+VWSSPEFCSCount);
	}
	
	
	//
	static class Script{
		//
		private boolean addCylindricalGS=false;
		
		private StringBuilder sb1=null;
		private StringBuilder sb2=null;
		
		private FileWriter fw=null;
		
		private CtlDescriptor ctl=null;
		
		
		//
		public Script(CtlDescriptor ctl,String path){
			this.ctl =ctl ;
			
			try{ fw=new FileWriter(path);}
			catch(IOException e){ e.printStackTrace(); System.exit(0);}
			
			sb1=new StringBuilder();
			sb1.append("'open "+ctl.getPath()+"'\n");
			sb1.append("'enable print "+path.replace(".gs",".gmf")+"'\n\n");
			
			sb2=new StringBuilder();
			sb2.append("\n\n'open "+ctl.getPath()+"'\n");
			sb2.append("'set gxout fwrite'\n");
			sb2.append("'set fwrite "+path.replace(".gs",".dat")+"'\n\n");
		}
		
		
		public void addDrawScript(Typhoon tr,int l,float[] E,float[] V,float[] S,float[] M_W){
			int xctr=ctl.getXLENum(tr.getLongitudes()[l]);
			int yctr=ctl.getYLENum(tr.getLatitudes()[l]);
			
			int radx=10;
			int rady= 8;
			
			int xstr=xctr-radx+1,xend=xctr+radx+2;
			int ystr=yctr-rady+1,yend=yctr+rady+2;
			
			float wholeU=tr.getZonalVelocity()[l];
			float wholeV=tr.getMeridionalVelocity()[l];
			
			float[] rw=tr.getStormRelativeWinds();
			float[] pr=tr.getPressures();
			
			MDate str=new MDate(tr.getTime(0));
			MDate tim=new MDate(tr.getTime(l));
			
			String ctlname=str.getYear()+String.format("%02d",str.getMonth())+tr.getName()+".ctl";
			
			if(addCylindricalGS){
				sb1.append("'open "+respath+"Station/"+ctlname+"'\n");
				sb1.append("'open "+respath+"Station/grid.ctl'\n");
			}
			sb1.append("'setvpage 2 2.4 2 1'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 200'\n");
			sb1.append("'set x "+xstr+" "+xend+"'\n");
			sb1.append("'set y "+ystr+" "+yend+"'\n");
			sb1.append("'set time "+tim.toGradsDate()+"'\n");
			if(addCylindricalGS) addCylindricalGS(l,E);
			else addPVTGS();
			sb1.append("'set gxout contour'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set ccolor 1'\n");
			sb1.append("'set arrscl 0.6 80'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 80'\n");
			sb1.append("'d u;v'\n");
			sb1.append("'drawmark 3 "+tr.getLongitudes()[l]+" "+tr.getLatitudes()[l]+" 0.2 '\n");
			sb1.append("'draw title "+(tr.getName()==null?"":tr.getName())+" "+
				tim.toGradsDate()+
				"\\P:"  +String.format("%.1f",pr[l])+
				" dP:" +String.format("%.1f",Typhoon.getChangesByCentralDiff(pr)[l])+
				" V:"  +String.format("%.1f",rw[l])+
				" dV:" +String.format("%.1f",Typhoon.getChangesByCentralDiff(rw)[l])+
				" M-W:"+String.format("%.1f",M_W[l])+"'\n\n"
			);
			sb1.append("'setvpage 2 2.4 1 1'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 850'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set arrscl 0.6 40'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 40'\n");
			sb1.append("'set rbcols auto'\n");
			sb1.append("'d u;v;mag(u,v)'\n");
			sb1.append("'drawmark 3 "+tr.getLongitudes()[l]+" "+tr.getLatitudes()[l]+" 0.2 '\n");
			sb1.append("'draw title" +
				" EFC:" +String.format("%.1f",E[l])+
				" VWS:" +String.format("%.1f",V[l])+
				" SST:" +String.format("%.1f",S[l])+"'\n\n"
			);
			sb1.append("'setvpage 2 2.4 2 2'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 200'\n");
			sb1.append("'set x "+xstr+" "+xend+"'\n");
			sb1.append("'set y "+ystr+" "+yend+"'\n");
			sb1.append("'set time "+tim.toGradsDate()+"'\n");
			if(addCylindricalGS) addCylindricalGS(l,E);
			else addPVTGS();
			sb1.append("'set gxout contour'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set ccolor 1'\n");
			sb1.append("'set arrscl 0.6 80'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 80'\n");
			sb1.append("'d u-"+wholeU+";v-"+wholeV+"'\n");
			sb1.append("'drawmark 3 "+tr.getLongitudes()[l]+" "+tr.getLatitudes()[l]+" 0.2 '\n");
			sb1.append("'draw title "+(tr.getName()==null?"":tr.getName())+" "+
				tim.toGradsDate()+
				"\\P:"  +String.format("%.1f",pr[l])+
				" dP:" +String.format("%.1f",Typhoon.getChangesByCentralDiff(pr)[l])+
				" V:"  +String.format("%.1f",rw[l])+
				" dV:" +String.format("%.1f",Typhoon.getChangesByCentralDiff(rw)[l])+
				" M-W:"+String.format("%.1f",M_W[l])+"'\n\n"
			);
			sb1.append("'setvpage 2 2.4 1 2'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 850'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set arrscl 0.6 40'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 40'\n");
			sb1.append("'set rbcols auto'\n");
			sb1.append("'d u-"+wholeU+";v-"+wholeV+";mag(u-"+wholeU+",v-"+wholeV+")'\n");
			sb1.append("'drawmark 3 "+tr.getLongitudes()[l]+" "+tr.getLatitudes()[l]+" 0.2 '\n");
			sb1.append("'draw title" +
				" EFC:" +String.format("%.1f",E[l])+
				" VWS:" +String.format("%.1f",V[l])+
				" SST:" +String.format("%.1f",S[l])+"'\n\n"
			);
			sb1.append("'print'\n");
			sb1.append("'c'\n");
			if(addCylindricalGS){
				sb1.append("'close 3'\n");
				sb1.append("'close 2'\n\n");
			}
			
			sb2.append("'set x "+xstr+" "+xend+"'\n");
			sb2.append("'set y "+ystr+" "+yend+"'\n");
			sb2.append("'set time "+tim.toGradsDate()+"'\n");
			sb2.append("'set lev 850'\n");
			sb2.append("'d u'\n");
			sb2.append("'set lev 200'\n");
			sb2.append("'d u'\n");
			sb2.append("'set lev 850'\n");
			sb2.append("'d v'\n");
			sb2.append("'set lev 200'\n");
			sb2.append("'d v'\n\n");
			sb2.append("'set lev 850'\n");
			sb2.append("'d u-"+wholeU+"'\n");
			sb2.append("'set lev 200'\n");
			sb2.append("'d u-"+wholeU+"'\n");
			sb2.append("'set lev 850'\n");
			sb2.append("'d v-"+wholeV+"'\n");
			sb2.append("'set lev 200'\n");
			sb2.append("'d v-"+wholeV+"'\n\n");
		}
		
		public void addCylindricalGS(int l,float[] E){
			boolean posEFC=E[l]>0;
			
			sb1.append("'set gxout shaded'\n");
			if(posEFC){
				sb1.append("'set rbcols 6 2 8 12 7 10 3 13 5 11 4 14 9'\n");
				sb1.append("'set black -200 200'\n");
				sb1.append("'set cint 200'\n");
				sb1.append("'set rbrange -2400 400'\n");
			}else{
				sb1.append("'set black -150 150'\n");
				sb1.append("'set cint 150'\n");
				sb1.append("'set rbrange -300 1800'\n");
			}
			sb1.append("'d oacres(u.3,flx.2)'\n");
			sb1.append("'cbarn 1 1 8.5 4.5'\n");
		}
		
		public void addPVTGS(){
			sb1.append("'set gxout shaded'\n");
			sb1.append("'set cint 2'\n");
			sb1.append("'set rbrange 212 232'\n");
			sb1.append("'d t'\n");
			sb1.append("'cbarn 1 1 8.5 4.5'\n");
			sb1.append("'set gxout contour'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set cint 1'\n");
			sb1.append("'d pv*1e6'\n");
		}
		
		public void flush(){
			sb1.append("'close 1'\n");
			sb1.append("'disable print'\n");
			sb1.append("'reinit'\n");
			
			sb2.append("'close 1'\n");
			sb2.append("'disable fwrite'\n");
			sb2.append("'reinit'\n");
			
			try{
				fw.write(sb1.toString());
				fw.write(sb2.toString());
				fw.close();
				
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
	}
}
