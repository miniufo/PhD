//
package Package;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;
import static Package.IntensityModel.*;


//
public class CompositeEFC{
	//
	private static int[] PEFCSWW=new int[3];
	private static int[] PEFCSNW=new int[3];
	private static int[] PEFCSNN=new int[3];
	private static int[] PEFCSNE=new int[3];
	private static int[] PEFCSEE=new int[3];
	private static int[] PEFCSSE=new int[3];
	private static int[] PEFCSSS=new int[3];
	private static int[] PEFCSSW=new int[3];
	
	private static int[] NEFCSWW=new int[3];
	private static int[] NEFCSNW=new int[3];
	private static int[] NEFCSNN=new int[3];
	private static int[] NEFCSNE=new int[3];
	private static int[] NEFCSEE=new int[3];
	private static int[] NEFCSSE=new int[3];
	private static int[] NEFCSSS=new int[3];
	private static int[] NEFCSSW=new int[3];
	
	private static final float threshold=10f;
	
	private static final boolean noDepress=false;	// wind > 17.2 m/s
	private static final boolean noLanding=false;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	private static final boolean forwardDf=true;	// using forward difference to get deltaP
	
	private static final StringBuilder sb=new StringBuilder();
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Statistics/"+dsets+"/Composite/";
	private static final String tranges="time=1Jan1987-31Dec2011";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		Script ComPEFCS_WW=new Script(ctl,respath+"ComPEFCS_WW.gs"); List<float[][][]> P_WW=new ArrayList<>();
		Script ComPEFCS_NW=new Script(ctl,respath+"ComPEFCS_NW.gs"); List<float[][][]> P_NW=new ArrayList<>();
		Script ComPEFCS_NN=new Script(ctl,respath+"ComPEFCS_NN.gs"); List<float[][][]> P_NN=new ArrayList<>();
		Script ComPEFCS_NE=new Script(ctl,respath+"ComPEFCS_NE.gs"); List<float[][][]> P_NE=new ArrayList<>();
		Script ComPEFCS_EE=new Script(ctl,respath+"ComPEFCS_EE.gs"); List<float[][][]> P_EE=new ArrayList<>();
		Script ComPEFCS_SE=new Script(ctl,respath+"ComPEFCS_SE.gs"); List<float[][][]> P_SE=new ArrayList<>();
		Script ComPEFCS_SS=new Script(ctl,respath+"ComPEFCS_SS.gs"); List<float[][][]> P_SS=new ArrayList<>();
		Script ComPEFCS_SW=new Script(ctl,respath+"ComPEFCS_SW.gs"); List<float[][][]> P_SW=new ArrayList<>();
		
		Script ComNEFCS_WW=new Script(ctl,respath+"ComNEFCS_WW.gs");
		Script ComNEFCS_NW=new Script(ctl,respath+"ComNEFCS_NW.gs");
		Script ComNEFCS_NN=new Script(ctl,respath+"ComNEFCS_NN.gs");
		Script ComNEFCS_NE=new Script(ctl,respath+"ComNEFCS_NE.gs");
		Script ComNEFCS_EE=new Script(ctl,respath+"ComNEFCS_EE.gs");
		Script ComNEFCS_SE=new Script(ctl,respath+"ComNEFCS_SE.gs");
		Script ComNEFCS_SS=new Script(ctl,respath+"ComNEFCS_SS.gs");
		Script ComNEFCS_SW=new Script(ctl,respath+"ComNEFCS_SW.gs");
		
		int cc=0,tycount=0,samples=0;
		for(Typhoon tr:ls){
			if(++cc%20==0) System.out.print(".");
			
			tycount++;
			samples+=tr.getTCount();
			
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,28,2,0.3f,-650,850));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=new Variable[3];
			Variable[] tmp=df.getVariables(new Range("",dd),false,"u","v");
			vars[0]=tmp[0];
			vars[1]=tmp[1];
			vars[2]=df.getVariables(new Range("z(1,1)",dd),false,"sst")[0];
			
			Variable[] shrs=dm.cVerticalWindShear(vars[0],vars[1]);
			
			Variable shrsum=dm.cRadialAverage(shrs[0],1,15).anomalizeX();	// 500 km
			Variable shrsvm=dm.cRadialAverage(shrs[1],1,15).anomalizeX();	// 500 km
			
			Variable uave=dm.cRadialAverage(vars[0],1,15).anomalizeX();		// 500 km
			Variable vave=dm.cRadialAverage(vars[1],1,15).anomalizeX();		// 500 km
			
			Variable relU=vars[0].copy();
			Variable relV=vars[1].copy();
			cRelativeVelocity(relU,relV,tr.getZonalVelocity(),tr.getMeridionalVelocity());
			Variable uthre=dm.cRadialAverage(relU,1,15).anomalizeX();		// 500 km
			Variable vthre=dm.cRadialAverage(relV,1,15).anomalizeX();		// 500 km
			
			Variable vwsm=shrsum.hypotenuse(shrsvm);
			Variable sstm=dm.cRadialAverage(vars[2],1,15).anomalizeX().minusEq(273.15f);
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getZonalVelocity(),tr.getMeridionalVelocity(),utvr[0],utvr[1]);
			
			Variable utm=utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable flux =utvr[0].multiply(utvr[1]);      // flux u'v'
			Variable lcefc=dm.cLocalREFC(utvr[0],utvr[1]); // localized REFC
			Variable iner =dm.cInertialStability(utvr[0],utvr[1]);	// local inertial stability
			Variable inerN=dm.cInertialStabilityNorm(utvr[0],utvr[1]);	// local inertial stability normalized by f2
			Variable absW =dm.cAbsoluteAngularVelocity(utvr[0]);	// local absolute angular velocity
			Variable eta  =dm.cAbsoluteVorticity(utvr[0],utvr[1]);	// local absolute vorticity
			Variable etam =dm.cMeanAbsoluteVorticity(utm);
			Variable efcsm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y,9,18);	// 300-600 km
			Variable isbm =dm.cMeanInertialStabilityByUT(utm).averageAlong(Dimension.Y,9,18);	// 300-600 km
			
			ct=new CoordinateTransformation(new SphericalSpatialModel(ctl2),csm);
			Variable lsmm=dm.cRadialAverage(ct.transToCylindricalInvariantly(lsm),1,6).anomalizeX();
			
			float[] ETAM= etam.getData()[1][0][0];
			float[] EFCS=efcsm.getData()[1][0][0];
			float[] ISBM= isbm.getData()[1][0][0];
			float[] VWSM= vwsm.getData()[0][0][0];
			float[] uppU= uave.getData()[1][0][0];
			float[] uppV= vave.getData()[1][0][0];
			float[] lowU= uave.getData()[0][0][0];
			float[] lowV= vave.getData()[0][0][0];
			float[] thrU=uthre.getData()[1][0][0];
			float[] thrV=vthre.getData()[1][0][0];
			float[] SSTM= sstm.getData()[0][0][0];
			float[] pote=cPotential(sstm.getData()[0][0][0],tr.getStormRelativeWinds());
			float[] degW=cAtanByWind(thrU,thrV);
			
			boolean[] wind=noDepress?greaterEqualThan(tr.getWinds(),17.2f):newBooleans(tr.getTCount());
			boolean[] lsmb=noLanding?lessThan(lsmm.getData()[0][0][0],1e-9f):newBooleans(tr.getTCount());
			boolean[] llon=withinWNP?greaterEqualThan(tr.getLongitudes(),100):newBooleans(tr.getTCount());
			boolean[] rlon=withinWNP?lessEqualThan(tr.getLongitudes(),190):newBooleans(tr.getTCount());
			
			boolean[] PEFCS=greaterEqualThan(efcsm.getData()[1][0][0],threshold);
			boolean[] NEFCS=lessEqualThan(efcsm.getData()[1][0][0],-threshold);
			
			boolean[] combPEFCS_WW=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW,   0,22.5f));
			boolean[] combPEFCS_NW=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW, -45,22.5f));
			boolean[] combPEFCS_NN=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW, -90,22.5f));
			boolean[] combPEFCS_NE=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW,-135,22.5f));
			boolean[] combPEFCS_EE=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW,-180,22.5f));
			boolean[] combPEFCS_SE=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW, 135,22.5f));
			boolean[] combPEFCS_SS=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW,  90,22.5f));
			boolean[] combPEFCS_SW=combination(wind,PEFCS,lsmb,llon,rlon,inRange(degW,  45,22.5f));
			
			boolean[] combNEFCS_WW=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW,   0,22.5f));
			boolean[] combNEFCS_NW=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW, -45,22.5f));
			boolean[] combNEFCS_NN=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW, -90,22.5f));
			boolean[] combNEFCS_NE=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW,-135,22.5f));
			boolean[] combNEFCS_EE=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW,-180,22.5f));
			boolean[] combNEFCS_SE=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW, 135,22.5f));
			boolean[] combNEFCS_SS=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW,  90,22.5f));
			boolean[] combNEFCS_SW=combination(wind,NEFCS,lsmb,llon,rlon,inRange(degW,  45,22.5f));
			
			float[] deltaP=forwardDf?
				Typhoon.getChangesByForwardDiff(tr.getPressures(),1):
				Typhoon.getChangesByCentralDiff(tr.getPressures());
			
			float[] deltaV=forwardDf?
				Typhoon.getChangesByForwardDiff(tr.getStormRelativeWinds(),1):
				Typhoon.getChangesByCentralDiff(tr.getStormRelativeWinds());
			
			for(int l=0;l<tr.getTCount();l++){
				if(EFCS[l]>20&&VWSM[l]<10&&deltaP[l]!=-9999f){
					String group=null;
					
					if(group==null&&combPEFCS_WW[l]) group="WW";
					if(group==null&&combPEFCS_NW[l]) group="NW";
					if(group==null&&combPEFCS_NN[l]) group="NN";
					if(group==null&&combPEFCS_NE[l]) group="NE";
					if(group==null&&combPEFCS_EE[l]) group="EE";
					if(group==null&&combPEFCS_SE[l]) group="SE";
					if(group==null&&combPEFCS_SS[l]) group="SS";
					if(group==null&&combPEFCS_SW[l]) group="SW";
					
					sb.append(String.format(
						"%10s  %5s  %s  %4s  EFC %7.2f    Eta %7.2f    VWS %7.2f    SST %7.2f    Pr %7.2f    dP %7.2f    dV %7.2f    M_W %7.2f    E/V %7.2f    Group %s\n",
						tr.getName(),tr.getID(),tr.getRecord(l),tr.getTypes()[l],EFCS[l],ETAM[l]*1e5f,VWSM[l],SSTM[l],tr.getPressures()[l],deltaP[l],deltaV[l],pote[l],EFCS[l]/VWSM[l],group
					));
				}
				
				// positive
				if(combPEFCS_WW[l]){
					ComPEFCS_WW.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_WW.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSWW[0]++;
					else if(dp==0) PEFCSWW[1]++;
					else PEFCSWW[2]++;
				}
				if(combPEFCS_NW[l]){
					ComPEFCS_NW.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_NW.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSNW[0]++;
					else if(dp==0) PEFCSNW[1]++;
					else PEFCSNW[2]++;
				}
				if(combPEFCS_NN[l]){
					ComPEFCS_NN.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_NN.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSNN[0]++;
					else if(dp==0) PEFCSNN[1]++;
					else PEFCSNN[2]++;
				}
				if(combPEFCS_NE[l]){
					ComPEFCS_NE.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_NE.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSNE[0]++;
					else if(dp==0) PEFCSNE[1]++;
					else PEFCSNE[2]++;
				}
				if(combPEFCS_EE[l]){
					ComPEFCS_EE.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_EE.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSEE[0]++;
					else if(dp==0) PEFCSEE[1]++;
					else PEFCSEE[2]++;
				}
				if(combPEFCS_SE[l]){
					ComPEFCS_SE.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_SE.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSSE[0]++;
					else if(dp==0) PEFCSSE[1]++;
					else PEFCSSE[2]++;
				}
				if(combPEFCS_SS[l]){
					ComPEFCS_SS.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_SS.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSSS[0]++;
					else if(dp==0) PEFCSSS[1]++;
					else PEFCSSS[2]++;
				}
				if(combPEFCS_SW[l]){
					ComPEFCS_SW.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					P_SW.add(getCylindricalData(l,relU,relV,flux,lcefc,iner,inerN,absW,eta));
					
					float dp=deltaP[l];
					if(dp<0) PEFCSSW[0]++;
					else if(dp==0) PEFCSSW[1]++;
					else PEFCSSW[2]++;
				}
				
				// negative
				if(combNEFCS_WW[l]){
					ComNEFCS_WW.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSWW[0]++;
					else if(dp==0) NEFCSWW[1]++;
					else NEFCSWW[2]++;
				}
				if(combNEFCS_NW[l]){
					ComNEFCS_NW.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSNW[0]++;
					else if(dp==0) NEFCSNW[1]++;
					else NEFCSNW[2]++;
				}
				if(combNEFCS_NN[l]){
					ComNEFCS_NN.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSNN[0]++;
					else if(dp==0) NEFCSNN[1]++;
					else NEFCSNN[2]++;
				}
				if(combNEFCS_NE[l]){
					ComNEFCS_NE.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSNE[0]++;
					else if(dp==0) NEFCSNE[1]++;
					else NEFCSNE[2]++;
				}
				if(combNEFCS_EE[l]){
					ComNEFCS_EE.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSEE[0]++;
					else if(dp==0) NEFCSEE[1]++;
					else NEFCSEE[2]++;
				}
				if(combNEFCS_SE[l]){
					ComNEFCS_SE.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSSE[0]++;
					else if(dp==0) NEFCSSE[1]++;
					else NEFCSSE[2]++;
				}
				if(combNEFCS_SS[l]){
					ComNEFCS_SS.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSSS[0]++;
					else if(dp==0) NEFCSSS[1]++;
					else NEFCSSS[2]++;
				}
				if(combNEFCS_SW[l]){
					ComNEFCS_SW.addDrawScript(tr,l,ETAM,ISBM,EFCS,VWSM,SSTM,pote,uppU,uppV,lowU,lowV,deltaP,deltaV);
					
					float dp=deltaP[l];
					if(dp<0) NEFCSSW[0]++;
					else if(dp==0) NEFCSSW[1]++;
					else NEFCSSW[2]++;
				}
			}
		}
		
		System.out.println(sb.toString());
		
		ComPEFCS_WW.flush();	ComPEFCS_NW.flush();	ComPEFCS_NN.flush();	ComPEFCS_NE.flush();
		ComPEFCS_EE.flush();	ComPEFCS_SE.flush();	ComPEFCS_SS.flush();	ComPEFCS_SW.flush();
		
		ComNEFCS_WW.flush();	ComNEFCS_NW.flush();	ComNEFCS_NN.flush();	ComNEFCS_NE.flush();
		ComNEFCS_EE.flush();	ComNEFCS_SE.flush();	ComNEFCS_SS.flush();	ComNEFCS_SW.flush();
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("forward diff for delP: "+forwardDf);
		System.out.println("\nthere are "+tycount+" TCs and "+samples+" 6hr samples");
		System.out.println("\nthreshold: "+threshold+" m s^-1 / day\n");
		
		System.out.println("PEFCS WW count: "+format(PEFCSWW)+ComPEFCS_WW.getMeanStatistics());
		System.out.println("PEFCS NW count: "+format(PEFCSNW)+ComPEFCS_NW.getMeanStatistics());
		System.out.println("PEFCS NN count: "+format(PEFCSNN)+ComPEFCS_NN.getMeanStatistics());
		System.out.println("PEFCS NE count: "+format(PEFCSNE)+ComPEFCS_NE.getMeanStatistics());
		System.out.println("PEFCS EE count: "+format(PEFCSEE)+ComPEFCS_EE.getMeanStatistics());
		System.out.println("PEFCS SE count: "+format(PEFCSSE)+ComPEFCS_SE.getMeanStatistics());
		System.out.println("PEFCS SS count: "+format(PEFCSSS)+ComPEFCS_SS.getMeanStatistics());
		System.out.println("PEFCS SW count: "+format(PEFCSSW)+ComPEFCS_SW.getMeanStatistics());
		System.out.println();
		System.out.println("NEFCS WW count: "+format(NEFCSWW)+ComNEFCS_WW.getMeanStatistics());
		System.out.println("NEFCS NW count: "+format(NEFCSNW)+ComNEFCS_NW.getMeanStatistics());
		System.out.println("NEFCS NN count: "+format(NEFCSNN)+ComNEFCS_NN.getMeanStatistics());
		System.out.println("NEFCS NE count: "+format(NEFCSNE)+ComNEFCS_NE.getMeanStatistics());
		System.out.println("NEFCS EE count: "+format(NEFCSEE)+ComNEFCS_EE.getMeanStatistics());
		System.out.println("NEFCS SE count: "+format(NEFCSSE)+ComNEFCS_SE.getMeanStatistics());
		System.out.println("NEFCS SS count: "+format(NEFCSSS)+ComNEFCS_SS.getMeanStatistics());
		System.out.println("NEFCS SW count: "+format(NEFCSSW)+ComNEFCS_SW.getMeanStatistics());
		
		writeCylindricalData(P_WW,respath+"CylindP_WW.dat");
		writeCylindricalData(P_NW,respath+"CylindP_NW.dat");
		writeCylindricalData(P_NN,respath+"CylindP_NN.dat");
		writeCylindricalData(P_NE,respath+"CylindP_NE.dat");
		writeCylindricalData(P_EE,respath+"CylindP_EE.dat");
		writeCylindricalData(P_SE,respath+"CylindP_SE.dat");
		writeCylindricalData(P_SS,respath+"CylindP_SS.dat");
		writeCylindricalData(P_SW,respath+"CylindP_SW.dat");
	}
	
	public static void cRelativeVelocity(Variable u,Variable v,float[] su,float[] sv){
		int t=u.getTCount();
		int z=u.getZCount();
		int y=u.getYCount();
		int x=u.getXCount();
		
		float[][][][] udata=u.getData();
		float[][][][] vdata=v.getData();
		
		for(int l=0;l<t;l++)
		for(int k=0;k<z;k++)
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			udata[k][j][i][l]-=su[l];
			vdata[k][j][i][l]-=sv[l];
		}
	}
	
	public static String format(int[] data){
		int sum=0;
		
		for(int d:data) sum+=d;
		
		return String.format("%3d [%3d(%5.1f), %3d(%5.1f), %3d(%5.1f)]",
			sum,
			data[0],sum==0?0:data[0]/(float)sum*100,
			data[1],sum==0?0:data[1]/(float)sum*100,
			data[2],sum==0?0:data[2]/(float)sum*100);
	}
	
	public static float[][][] getCylindricalData
	(int l,Variable relU,Variable relV,Variable flux,Variable lcefc,Variable iner,Variable inerN,Variable absW,Variable eta){
		int y=relU.getYCount(),x=relV.getXCount();
		
		float[][][] re=new float[8][y][x];
		float[][][] ud=relU.getData()[1];
		float[][][] vd=relV.getData()[1];
		float[][][] fd=flux.getData()[1];
		float[][][] ld=lcefc.getData()[1];
		float[][][] id=iner.getData()[1];
		float[][][] nd=inerN.getData()[1];
		float[][][] ad=absW.getData()[1];
		float[][][] ed=eta.getData()[1];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			re[0][j][i]=ud[j][i][l];
			re[1][j][i]=vd[j][i][l];
			re[2][j][i]=fd[j][i][l];
			re[3][j][i]=ld[j][i][l];
			re[4][j][i]=id[j][i][l];
			re[5][j][i]=nd[j][i][l];
			re[6][j][i]=ad[j][i][l];
			re[7][j][i]=ed[j][i][l];
		}
		
		return re;
	}
	
	public static void writeCylindricalData(List<float[][][]> data,String fname){
		int t=data.size(),y=data.get(0)[0].length,x=data.get(0)[0][0].length;
		
		StringBuilder sb=new StringBuilder();
		
		sb.append(
			"dpath ^Haima.ctl\n"+
			"title haima (200421) CMA Data\n"+
			"xdef  "+x+"\n"+
			"ydef  "+y+" 0.3\n"+
			"zdef   2 850  -650\n"+
			"tdef  "+t+" 00z11sep2004 6hr\n"+
			"coords\n"			
		);
		
		for(int l=0;l<t;l++) sb.append("130 20 1000 15\n");
		
		sb.append("endcoords\n");
		
		Variable relU=new Variable("relU",new Range(t,1,y,x));
		Variable relV=new Variable("relV",new Range(t,1,y,x));
		Variable flux=new Variable("flux",new Range(t,1,y,x));
		Variable lefc=new Variable("lefc",new Range(t,1,y,x));
		Variable iner=new Variable("iner",new Range(t,1,y,x));
		Variable ineN=new Variable("inerN",new Range(t,1,y,x));
		Variable absW=new Variable("absW",new Range(t,1,y,x));
		Variable eta =new Variable("eta" ,new Range(t,1,y,x));
		
		relU.setUndef(-9.99e+08f); relU.setCommentAndUnit("relative U"); 
		relV.setUndef(-9.99e+08f); relV.setCommentAndUnit("relative V"); 
		flux.setUndef(-9.99e+08f); flux.setCommentAndUnit("flux of u'v'"); 
		lefc.setUndef(-9.99e+08f); lefc.setCommentAndUnit("local REFC");
		iner.setUndef(-9.99e+08f); iner.setCommentAndUnit("local inerital stability"); 
		ineN.setUndef(-9.99e+08f); ineN.setCommentAndUnit("normalized local inerital stability"); 
		absW.setUndef(-9.99e+08f); absW.setCommentAndUnit("absolute angular velocity"); 
		 eta.setUndef(-9.99e+08f);  eta.setCommentAndUnit("absolute vorticity"); 
		
		float[][][] udata=relU.getData()[0];
		float[][][] vdata=relV.getData()[0];
		float[][][] fdata=flux.getData()[0];
		float[][][] ldata=lefc.getData()[0];
		float[][][] idata=iner.getData()[0];
		float[][][] ndata=ineN.getData()[0];
		float[][][] adata=absW.getData()[0];
		float[][][] edata= eta.getData()[0];
		
		for(int l=0;l<t;l++)
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			udata[j][i][l]=data.get(l)[0][j][i];
			vdata[j][i][l]=data.get(l)[1][j][i];
			fdata[j][i][l]=data.get(l)[2][j][i];
			ldata[j][i][l]=data.get(l)[3][j][i];
			idata[j][i][l]=data.get(l)[4][j][i];
			ndata[j][i][l]=data.get(l)[5][j][i];
			adata[j][i][l]=data.get(l)[6][j][i];
			edata[j][i][l]=data.get(l)[7][j][i];
		}
		
		DiagnosisFactory df=DiagnosisFactory.parseContent(sb.toString());
		
		DataDescriptor dd=df.getDataDescriptor();
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,fname);
		dw.writeData(dd,relU,relV,flux,lefc,iner,ineN,absW,eta); dw.closeFile();
	}
	
	
	//
	static class Script{
		//
		private int count=0;
		private int countSST=0;
		
		private float aveSST=0;
		private float aveVWS=0;
		private float aveEFC=0;
		private float aveISB=0;
		private float aveETA=0;
		
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
		
		 
		public void addDrawScript(Typhoon tr,int l,float[] eta,float[] I,float[] E,float[] V,float[] S,float[] M_W,
		float[] uppU,float[] uppV,float[] lowU,float[] lowV,float[] dP,float[] dV){
			count++;
			
			if(S[l]>0){ aveSST+=S[l]; countSST++;}
			aveVWS+=V[l];
			aveEFC+=E[l];
			aveISB+=I[l]*1e9f;
			aveETA+=eta[l]*1e5f;
			
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
			
			MDate md=new MDate(tr.getTime(l));
			
			sb1.append("'setvpage 2 2.4 2 1'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 200'\n");
			sb1.append("'set x "+xstr+" "+xend+"'\n");
			sb1.append("'set y "+ystr+" "+yend+"'\n");
			sb1.append("'set time "+md.toGradsDate()+"'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set arrscl 0.6 80'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 80'\n");
			sb1.append("'d u;v;mag(u,v)'\n");
			sb1.append("'drawmark 3 "+tr.getLongitudes()[l]+" "+tr.getLatitudes()[l]+" 0.2 '\n");
			sb1.append("'draw title "+(tr.getName()==null?"":tr.getName())+" "+
				md.toGradsDate()+
				"\\P:"  +String.format("%.1f",pr[l])+
				" dP:" +String.format("%.1f",dP[l])+
				" V:"  +String.format("%.1f",rw[l])+
				" dV:" +String.format("%.1f",dV[l])+
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
				" SST:" +String.format("%.1f",S[l])+
				"\\aveU:"+String.format("%.1f %.1f",uppU[l],lowU[l])+
				" aveV:"+String.format("%.1f %.1f",uppV[l],lowV[l])+"'\n\n"
			);
			sb1.append("'setvpage 2 2.4 2 2'\n");
			sb1.append("'setlopts 10 0.24 5 5'\n");
			sb1.append("'set lev 200'\n");
			sb1.append("'set x "+xstr+" "+xend+"'\n");
			sb1.append("'set y "+ystr+" "+yend+"'\n");
			sb1.append("'set time "+md.toGradsDate()+"'\n");
			sb1.append("'set cthick 10'\n");
			sb1.append("'set arrscl 0.6 80'\n");
			sb1.append("'set arrowhead -0.35'\n");
			sb1.append("'set rbrange 0 80'\n");
			sb1.append("'d u-"+wholeU+";v-"+wholeV+";mag(u-"+wholeU+",v-"+wholeV+")'\n");
			sb1.append("'drawmark 3 "+tr.getLongitudes()[l]+" "+tr.getLatitudes()[l]+" 0.2 '\n");
			sb1.append("'draw title "+(tr.getName()==null?"":tr.getName())+" "+
				md.toGradsDate()+
				"\\P:"  +String.format("%.1f",pr[l])+
				" dP:" +String.format("%.1f",dP[l])+
				" V:"  +String.format("%.1f",rw[l])+
				" dV:" +String.format("%.1f",dV[l])+
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
				" SST:" +String.format("%.1f",S[l])+
				"\\aveU:"+String.format("%.1f %.1f",uppU[l],lowU[l])+
				" aveV:"+String.format("%.1f %.1f",uppV[l],lowV[l])+"'\n\n"
			);
			sb1.append("'print'\n");
			sb1.append("'c'\n");
			
			sb2.append("'set x "+xstr+" "+xend+"'\n");
			sb2.append("'set y "+ystr+" "+yend+"'\n");
			sb2.append("'set time "+md.toGradsDate()+"'\n");
			sb2.append("'set lev 200'\n");
			sb2.append("'d pv'\n");
			sb2.append("'d t'\n");
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
		
		public String getMeanStatistics(){
			return String.format(
				"; mean SST (%6.2f), VWS (%6.2f), EFC (%6.2f), ISB (%6.2f), ETA (%6.2f); sample: %3d",
				countSST==0?0:aveSST/countSST,
				count==0?0:aveVWS/count,
				count==0?0:aveEFC/count,
				count==0?0:aveISB/count,
				count==0?0:aveETA/count,
				count
			);
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
