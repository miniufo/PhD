//
package Package;

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
import miniufo.diagnosis.Variable.Dimension;
import miniufo.lagrangian.Typhoon;


//
public class ScriptEFC{
	//
	private static int PEFCSCount=0;
	private static int NEFCSCount=0;
	private static int PEFCLCount=0;
	private static int NEFCLCount=0;
	
	private static final float threshold=60f;
	
	private static final boolean noDepress=false;	// wind > 17.2 m/s
	private static final boolean noLanding=false;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	private static final boolean forwardDf=true;	// using forward difference to get deltaP
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Statistics/"+dsets+"/";
	private static final String tranges="time=1Jan1987-31Dec2011";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		Script scrpPEFCS=new Script(ctl,respath+"ScriptPEFCS.gs");
		Script scrpNEFCS=new Script(ctl,respath+"ScriptNEFCS.gs");
		Script scrpPEFCL=new Script(ctl,respath+"ScriptPEFCL.gs");
		Script scrpNEFCL=new Script(ctl,respath+"ScriptNEFCL.gs");
		
		int cc=0;
		for(Typhoon tr:ls){
			if(++cc%20==0) System.out.print(".");
			
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,25,2,0.3f,-650,850));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v","sst");
			Variable[] shrs=dm.cVerticalWindShear(vars[0],vars[1]);
			
			Variable shrsum=dm.cRadialAverage(shrs[0],1,15).anomalizeX();
			Variable shrsvm=dm.cRadialAverage(shrs[1],1,15).anomalizeX();
			
			Variable uave=dm.cRadialAverage(vars[0],1,15).anomalizeX();
			Variable vave=dm.cRadialAverage(vars[1],1,15).anomalizeX();
			
			Variable vwsm=shrsum.hypotenuse(shrsvm);
			Variable sstm=dm.cRadialAverage(vars[2],1,15).anomalizeX().minusEq(273.15f);
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
			Variable utm=utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efclm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),15,24);	// 500-800 km
			Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]), 9,18);	// 300-600 km
			Variable isbym=dm.cMeanInertialStabilityByUT(utm).averageAlong(Dimension.Y,9,18);
			Variable etam=dm.cMeanAbsoluteVorticity(utm).averageAlong(Dimension.Y,9,18);
			
			ct=new CoordinateTransformation(new SphericalSpatialModel(ctl2),csm);
			Variable lsmm=dm.cRadialAverage(ct.transToCylindricalInvariantly(lsm),1,6).anomalizeX();
			
			boolean[] wind=noDepress?	// wind >= 17.2 m/s
				IntensityModel.greaterEqualThan(tr.getWinds(),17.2f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] lsmb=noLanding?	// no land within 200 km
				IntensityModel.lessThan(lsmm.getData()[0][0][0],1e-9f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] llon=withinWNP?	// lons >= 100E
				IntensityModel.greaterEqualThan(tr.getXPositions(),100):IntensityModel.newBooleans(tr.getTCount());
			boolean[] rlon=withinWNP?	// lons <= 190E
				IntensityModel.lessEqualThan(tr.getXPositions(),190):IntensityModel.newBooleans(tr.getTCount());
			
			boolean[] hghPEFCS=IntensityModel.greaterEqualThan(efcsm.getData()[1][0][0],threshold+40);
			boolean[] hghNEFCS=IntensityModel.lessEqualThan(efcsm.getData()[1][0][0],-threshold);
			boolean[] hghPEFCL=IntensityModel.greaterEqualThan(efclm.getData()[1][0][0],threshold+40);
			boolean[] hghNEFCL=IntensityModel.lessEqualThan(efclm.getData()[1][0][0],-threshold);
			
			boolean[] combPEFCS=IntensityModel.combination(wind,hghPEFCS,lsmb,llon,rlon);
			boolean[] combNEFCS=IntensityModel.combination(wind,hghNEFCS,lsmb,llon,rlon);
			boolean[] combPEFCL=IntensityModel.combination(wind,hghPEFCL,lsmb,llon,rlon);
			boolean[] combNEFCL=IntensityModel.combination(wind,hghNEFCL,lsmb,llon,rlon);
			
			float[] EFCS=efcsm.getData()[1][0][0];
			float[] EFCL=efclm.getData()[1][0][0];
			float[] VWSM= vwsm.getData()[0][0][0];
			float[] uppU= uave.getData()[1][0][0];
			float[] uppV= vave.getData()[1][0][0];
			float[] lowU= uave.getData()[0][0][0];
			float[] lowV= vave.getData()[0][0][0];
			float[] SSTM= sstm.getData()[0][0][0];
			float[] isbm=isbym.getData()[1][0][0];
			float[] meta= etam.getData()[1][0][0];
			float[] pote=IntensityModel.cPotential(sstm.getData()[0][0][0],tr.getStormRelativeWinds());
			
			for(int l=0;l<tr.getTCount();l++){
				if(combPEFCS[l]){ scrpPEFCS.addDrawScript(tr,l,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,isbm,meta);PEFCSCount++;}
				if(combNEFCS[l]){ scrpNEFCS.addDrawScript(tr,l,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,isbm,meta);NEFCSCount++;}
				if(combPEFCL[l]){ scrpPEFCL.addDrawScript(tr,l,EFCL,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,isbm,meta);PEFCLCount++;}
				if(combNEFCL[l]){ scrpNEFCL.addDrawScript(tr,l,EFCL,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,isbm,meta);NEFCLCount++;}
			}
		}
		
		scrpPEFCS.flush();	scrpNEFCS.flush();
		scrpPEFCL.flush();	scrpNEFCL.flush();
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("\nthreshold: "+threshold+" m s^-1 / day");
		System.out.println("effective PEFCS count:"+PEFCSCount);
		System.out.println("effective NEFCS count:"+NEFCSCount);
		System.out.println("effective PEFCL count:"+PEFCLCount);
		System.out.println("effective NEFCL count:"+NEFCLCount);
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
		
		
		public void addDrawScript(Typhoon tr,int l,float[] E,float[] V,float[] S,float[] M_W,
		float[] uppU,float[] uppV,float[] lowU,float[] lowV,float[] isbm,float[] etam){
			int xctr=ctl.getXLENum(tr.getXPositions()[l]);
			int yctr=ctl.getYLENum(tr.getYPositions()[l]);
			
			int radx=10;
			int rady= 8;
			
			int xstr=xctr-radx+1,xend=xctr+radx+2;
			int ystr=yctr-rady+1,yend=yctr+rady+2;
			
			float wholeU=tr.getUVel()[l];
			float wholeV=tr.getVVel()[l];
			
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
			//if(addCylindricalGS) addCylindricalGS(l,E);
			//else addPVTGS();
			sb1.append("'set gxout shade2'\n");
			sb1.append("'set cmin 40'\n");
			sb1.append("'set rbrange 40 200'\n");
			sb1.append("'set cint 40'\n");
			sb1.append("'d mag(u,v)'\n");
			sb1.append("'set gxout contour'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set ccolor 1'\n");
			sb1.append("'set arrscl 0.6 80'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 80'\n");
			sb1.append("'d u;v'\n");
			sb1.append("'drawmark 3 "+tr.getXPositions()[l]+" "+tr.getYPositions()[l]+" 0.2 '\n");
			sb1.append("'draw title "+(tr.getName()==null?"":tr.getName())+" "+
				tim.toGradsDate()+
				"\\P:"  +String.format("%.1f",pr[l])+
				" dP:" +String.format("%.1f",forwardDf?Typhoon.getChangesByForwardDiff(pr,1)[l]:Typhoon.getChangesByCentralDiff(pr)[l])+
				" V:"  +String.format("%.1f",rw[l])+
				" dV:" +String.format("%.1f",forwardDf?Typhoon.getChangesByForwardDiff(rw,1)[l]:Typhoon.getChangesByCentralDiff(rw)[l])+
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
			sb1.append("'drawmark 3 "+tr.getXPositions()[l]+" "+tr.getYPositions()[l]+" 0.2 '\n");
			sb1.append("'draw title" +
				" EFC:" +String.format("%.1f",E[l])+
				" ISB:" +String.format("%.3f",isbm[l]*1e9)+
				" ETA:" +String.format("%.3f",etam[l]*1e5)+
				" VWS:" +String.format("%.1f",V[l])+
				" SST:" +String.format("%.1f",S[l])+
				"\\aveU:"+String.format("%.1f %.1f",uppU[l],lowU[l])+
				" aveV:"+String.format("%.1f %.1f",uppV[l],lowV[l])+"'\n\n"
			);
			sb1.append("'setvpage 2 2.4 2 2'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 200'\n");
			sb1.append("'set x "+xstr+" "+xend+"'\n");
			sb1.append("'set y "+ystr+" "+yend+"'\n");
			sb1.append("'set time "+tim.toGradsDate()+"'\n");
			//if(addCylindricalGS) addCylindricalGS(l,E);
			//else addPVTGS();
			sb1.append("'set gxout shade2'\n");
			sb1.append("'set cmin 40'\n");
			sb1.append("'set rbrange 40 200'\n");
			sb1.append("'set cint 40'\n");
			sb1.append("'d mag(u-"+wholeU+",v-"+wholeV+")'\n");
			sb1.append("'set gxout contour'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set ccolor 1'\n");
			sb1.append("'set arrscl 0.6 80'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 80'\n");
			sb1.append("'d u-"+wholeU+";v-"+wholeV+"'\n");
			sb1.append("'drawmark 3 "+tr.getXPositions()[l]+" "+tr.getYPositions()[l]+" 0.2 '\n");
			sb1.append("'draw title "+(tr.getName()==null?"":tr.getName())+" "+
				tim.toGradsDate()+
				"\\P:"  +String.format("%.1f",pr[l])+
				" dP:" +String.format("%.1f",forwardDf?Typhoon.getChangesByForwardDiff(pr,1)[l]:Typhoon.getChangesByCentralDiff(pr)[l])+
				" V:"  +String.format("%.1f",rw[l])+
				" dV:" +String.format("%.1f",forwardDf?Typhoon.getChangesByForwardDiff(rw,1)[l]:Typhoon.getChangesByCentralDiff(rw)[l])+
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
			sb1.append("'drawmark 3 "+tr.getXPositions()[l]+" "+tr.getYPositions()[l]+" 0.2 '\n");
			sb1.append("'draw title" +
				" EFC:" +String.format("%.1f",E[l])+
				" ISB:" +String.format("%.3f",isbm[l]*1e9)+
				" ETA:" +String.format("%.3f",etam[l]*1e5)+
				" VWS:" +String.format("%.1f",V[l])+
				" SST:" +String.format("%.1f",S[l])+
				"\\aveU:"+String.format("%.1f %.1f",uppU[l],lowU[l])+
				" aveV:"+String.format("%.1f %.1f",uppV[l],lowV[l])+"'\n\n"
			);
			sb1.append("'print'\n");
			sb1.append("'c'\n");
			//if(addCylindricalGS){
			//	sb1.append("'close 3'\n");
			//	sb1.append("'close 2'\n\n");
			//}
			
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
