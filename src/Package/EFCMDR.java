package Package;

import java.util.List;

import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.lagrangian.Typhoon;

//
public class EFCMDR{
	//
	private static int ystr=1987;
	private static int yend=2011;
	
	private static int effectiveCount=0;
	private static int effectiveTC   =0;
	
	private static final boolean noDepress=false;	// wind > 17.2 m/s
	private static final boolean noLanding=false;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String tranges="time=1Jan"+ystr+"-31Dec"+yend;
	
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,25,2,0.3f,-650,850));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v");
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			//Variable efclm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y,15,24);	// 500-800 km
			Variable efcsm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y, 9,18);	// 300-600 km
			
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
			
			boolean[] EFCS=IntensityModel.ABSgreaterEqualThan(efcsm.getData()[1][0][0],10);// |EFCS| >= 10
			//boolean[] EFCL=IntensityModel.ABSgreaterEqualThan(efclm.getData()[1][0][0],10);// |EFCL| >= 10
			
			boolean[] validAll =IntensityModel.combination(wind,lsmb,llon,rlon);
			
			float[] lats=tr.getYPositions();
			
			int mi=getPeakIntensityTNum(tr.getPressures());
			
			//if(mi==1){ System.out.println(tr);System.exit(0);}
			
			for(int l=0,L=tr.getTCount();l<L;l++)
			if(lats[l]<20&&EFCS[l]) System.out.println(String.format("TYPE:%s     EFC:%8.2f     l: %3d     mi: %3d     L: %3d",
				tr.getTypes()[l],efcsm.getData()[1][0][0][l],l,mi,L
			));
			
			int cc=IntensityModel.getValidCount(validAll);
			effectiveCount+=cc;
			if(cc>0) effectiveTC++;
			//if(effectiveTC%20==0) System.out.print(".");
		}
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
	}
	
	static int getPeakIntensityTNum(float[] prs){
		float tmp=Float.MAX_VALUE;
		
		for(int l=0,L=prs.length;l<L;l++) if(tmp>prs[l]) tmp=prs[l];
		
		for(int l=0,L=prs.length;l<L;l++) if(prs[l]==tmp) return l+1;
		
		throw new IllegalArgumentException("no index returned");
	}
}
