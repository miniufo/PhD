//
package Package;

import java.util.ArrayList;
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.basic.ArrayUtil;
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
import miniufo.lagrangian.Typhoon.TYPE;
import miniufo.statistics.FilterModel;
import miniufo.statistics.StatisticsUtil;

//
public class EFCICPercentage{
	//
	private static int effectiveCount1=0;
	private static int effectiveTC1   =0;
	private static int effectiveCount2=0;
	private static int effectiveTC2   =0;
	
	private static float minLat= 90;
	private static float minLon=360;
	private static float maxLat=-90;
	private static float maxLon=  0;
	
	private static final boolean noDepress=true;	// wind > 17.2 m/s
	private static final boolean noLanding=true;	// no land within 100 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	private static final boolean hasPotent=true;
	
	private static final DataSets dsets=DataSets.JMA;
	//private static final String respath="d:/Data/PhD/Statistics/"+dsets+"/";
	private static final String tranges="time=1Jan1987-31Dec2011";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		List<InteractionEvent> res1=new ArrayList<InteractionEvent>();
		List<InteractionEvent> res2=new ArrayList<InteractionEvent>();
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,19,2,0.3f,-650,850));
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
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v","sst");
			Variable[] shrs=dm.cVerticalWindShear(vars[0],vars[1]);
			
			Variable shrsum=dm.cRadialAverage(shrs[0],1,15).anomalizeX();
			Variable shrsvm=dm.cRadialAverage(shrs[1],1,15).anomalizeX();
			
			Variable vwsm=shrsum.hypotenuse(shrsvm);
			Variable sstm=dm.cRadialAverage(vars[2],1,15).anomalizeX().minusEq(273.15f);
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			//Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),15,27);	// 500-900 km
			//Variable ffctm=dm.cRadialAverage(dm.cFFIndex(gaz,utvr[1]) ,15,27);	// 500-900 km
			
			/** 300-600 km */
			Variable efcsm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]), 9,18);	// 300-600 km
			
			ct=new CoordinateTransformation(new SphericalSpatialModel(ctl2),csm);
			Variable lsmm=dm.cRadialAverage(ct.transToCylindricalInvariantly(lsm),1,3).anomalizeX();
			
			boolean[] wind=noDepress?	// wind >= 17.2 m/s
				IntensityModel.greaterEqualThan(tr.getWinds(),17.2f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] lsmb=noLanding?	// no land within 200 km
				IntensityModel.lessThan(lsmm.getData()[0][0][0],1e-9f):IntensityModel.newBooleans(tr.getTCount());
			boolean[] llon=withinWNP?	// lons >= 100E
				IntensityModel.greaterEqualThan(tr.getXPositions(),100):IntensityModel.newBooleans(tr.getTCount());
			boolean[] rlon=withinWNP?	// lons <= 190E
				IntensityModel.lessEqualThan(tr.getXPositions(),190):IntensityModel.newBooleans(tr.getTCount());
			
			float[] EFC=efcsm.getData()[1][0][0];
			float[] VWS= vwsm.getData()[0][0][0];
			float[] SST= sstm.getData()[0][0][0];
			float[] POT=IntensityModel.cPotential(SST,tr.getWinds());
			
			//boolean[] hEFC=IntensityModel.lessEqualThan(EFC,-9f);
			//boolean[] hFFC=IntensityModel.lessEqualThan(FFC,-3f);
			boolean[] hEFC=IntensityModel.greaterEqualThan(EFC,10f);
			boolean[] hPOT=IntensityModel.greaterEqualThan(POT,30);
			
			boolean[] landing =IntensityModel.greaterEqualThan(lsmm.getData()[0][0][0],1e-9f);
			boolean[] validEFC=hasPotent?
				IntensityModel.combination(wind,lsmb,llon,rlon,hEFC,hPOT):
				IntensityModel.combination(wind,lsmb,llon,rlon,hEFC);
			boolean[] validPCH=hasPotent?
				IntensityModel.combination(wind,lsmb,llon,rlon,hPOT):
				IntensityModel.combination(wind,lsmb,llon,rlon);
			
			/**
			System.out.println("Times:           Lats     lons   Press    EFC   VWS   SST   POT   hEFC   hPCH Landing   hPOT");
			for(int l=0;l<tr.getCount();l++)
			System.out.println(String.format("%10d, %4.1f��N, %5.1f��E, %5.0f, %5.1f, %4.1f, %4.1f, %4.1f, %5s, %5s,  %5s, %5s",
				tr.getDates()[l].getLongTime(),tr.getLats()[l],tr.getLons()[l],tr.getPressures()[l],
				EFC[l],VWS[l],SST[l],POT[l],validEFC[l],validPCH[l],landing[l],hPOT[l]
			));*/
			
			InteractionEvent[] events1=getEvents(tr,validEFC,validPCH,landing,hPOT,EFC,SST,VWS,POT);
			
			if(events1!=null)
			for(InteractionEvent e:events1){ res1.add(e);} //System.out.println(e);}
			
			int cc=0;
			if(events1!=null) cc=events1.length;
			effectiveCount1+=cc;
			if(cc>0) effectiveTC1++;
		}
		
		int eventC=res1.size();
		
		Event[] events=new Event[eventC];
		
		for(int i=0;i<eventC;i++){
			InteractionEvent ie=res1.get(i);
			
			System.out.println(ie);
			
			events[i]=new Event(ie);
		}
		
		List<Event>[] re=groupEvents(events);
		
		compositeEvents(re[0]);
		compositeEvents(re[1]);
		
		System.exit(0);
		
		
		
		
		
		
		
		
		printResult(effectiveTC1,effectiveCount1);
		
		int ic_6h =NHourIntensityChangeCount(res1,1 );
		int ic_12h=NHourIntensityChangeCount(res1,2 );
		int ic_18h=NHourIntensityChangeCount(res1,3 );
		int ic_24h=NHourIntensityChangeCount(res1,4 );
		int ic_30h=NHourIntensityChangeCount(res1,5 );
		int ic_36h=NHourIntensityChangeCount(res1,6 );
		int ic_42h=NHourIntensityChangeCount(res1,7 );
		int ic_48h=NHourIntensityChangeCount(res1,8 );
		int ic_54h=NHourIntensityChangeCount(res1,9 );
		int ic_60h=NHourIntensityChangeCount(res1,10);
		int ic_66h=NHourIntensityChangeCount(res1,11);
		int ic_72h=NHourIntensityChangeCount(res1,12);
		
		System.out.println("\ncounts: "+
			ic_6h +"  "+ic_12h+"  "+ic_18h+"  "+ic_24h+"  "+
			ic_30h+"  "+ic_36h+"  "+ic_42h+"  "+ic_48h+"  "+
			ic_54h+"  "+ic_60h+"  "+ic_66h+"  "+ic_72h+"\n"
		);
		
		float[] EFC6M =new float[ic_6h ]; float[] pc_6h =new float[ic_6h ];
		float[] EFC12M=new float[ic_12h]; float[] pc_12h=new float[ic_12h];
		float[] EFC18M=new float[ic_18h]; float[] pc_18h=new float[ic_18h];
		float[] EFC24M=new float[ic_24h]; float[] pc_24h=new float[ic_24h];
		float[] EFC30M=new float[ic_30h]; float[] pc_30h=new float[ic_30h];
		float[] EFC36M=new float[ic_36h]; float[] pc_36h=new float[ic_36h];
		float[] EFC42M=new float[ic_42h]; float[] pc_42h=new float[ic_42h];
		float[] EFC48M=new float[ic_48h]; float[] pc_48h=new float[ic_48h];
		float[] EFC54M=new float[ic_54h]; float[] pc_54h=new float[ic_54h];
		float[] EFC60M=new float[ic_60h]; float[] pc_60h=new float[ic_60h];
		float[] EFC66M=new float[ic_66h]; float[] pc_66h=new float[ic_66h];
		float[] EFC72M=new float[ic_72h]; float[] pc_72h=new float[ic_72h];
		
		int pt_6h =0,pt_12h=0,pt_18h=0,pt_24h=0;
		int pt_30h=0,pt_36h=0,pt_42h=0,pt_48h=0;
		int pt_54h=0,pt_60h=0,pt_66h=0,pt_72h=0;
		for(InteractionEvent e:res1){
			if(e.hasIC(1 )){ EFC6M [pt_6h ]=e.cEFCMean(); pc_6h [pt_6h ]=e.IC(1 ); pt_6h ++;}
			if(e.hasIC(2 )){ EFC12M[pt_12h]=e.cEFCMean(); pc_12h[pt_12h]=e.IC(2 ); pt_12h++;}
			if(e.hasIC(3 )){ EFC18M[pt_18h]=e.cEFCMean(); pc_18h[pt_18h]=e.IC(3 ); pt_18h++;}
			if(e.hasIC(4 )){ EFC24M[pt_24h]=e.cEFCMean(); pc_24h[pt_24h]=e.IC(4 ); pt_24h++;}
			if(e.hasIC(5 )){ EFC30M[pt_30h]=e.cEFCMean(); pc_30h[pt_30h]=e.IC(5 ); pt_30h++;}
			if(e.hasIC(6 )){ EFC36M[pt_36h]=e.cEFCMean(); pc_36h[pt_36h]=e.IC(6 ); pt_36h++;}
			if(e.hasIC(7 )){ EFC42M[pt_42h]=e.cEFCMean(); pc_42h[pt_42h]=e.IC(7 ); pt_42h++;}
			if(e.hasIC(8 )){ EFC48M[pt_48h]=e.cEFCMean(); pc_48h[pt_48h]=e.IC(8 ); pt_48h++;}
			if(e.hasIC(9 )){ EFC54M[pt_54h]=e.cEFCMean(); pc_54h[pt_54h]=e.IC(9 ); pt_54h++;}
			if(e.hasIC(10)){ EFC60M[pt_60h]=e.cEFCMean(); pc_60h[pt_60h]=e.IC(10); pt_60h++;}
			if(e.hasIC(11)){ EFC66M[pt_66h]=e.cEFCMean(); pc_66h[pt_66h]=e.IC(11); pt_66h++;}
			if(e.hasIC(12)){ EFC72M[pt_72h]=e.cEFCMean(); pc_72h[pt_72h]=e.IC(12); pt_72h++;}
		}
		
		int[] tmp=null;	int[] tmp2=null;	float de=0;
		tmp=countIC(pc_6h);		tmp2=countICAmount(pc_6h,1);	de=ic_6h /100f;
		System.out.println(String.format(" 6h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_6h ,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC6M)));
		tmp=countIC(pc_12h);	tmp2=countICAmount(pc_12h,2);	de=ic_12h/100f;
		System.out.println(String.format("12h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_12h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC12M)));
		tmp=countIC(pc_18h);	tmp2=countICAmount(pc_18h,3);	de=ic_18h/100f;
		System.out.println(String.format("18h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_18h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC18M)));
		tmp=countIC(pc_24h);	tmp2=countICAmount(pc_24h,4);	de=ic_24h/100f;
		System.out.println(String.format("24h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_24h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC24M)));
		tmp=countIC(pc_30h);	tmp2=countICAmount(pc_30h,5);	de=ic_30h/100f;
		System.out.println(String.format("30h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_30h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC30M)));
		tmp=countIC(pc_36h);	tmp2=countICAmount(pc_36h,6);	de=ic_36h/100f;
		System.out.println(String.format("36h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_36h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC36M)));
		tmp=countIC(pc_42h);	tmp2=countICAmount(pc_42h,7);	de=ic_42h/100f;
		System.out.println(String.format("42h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_42h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC42M)));
		tmp=countIC(pc_48h);	tmp2=countICAmount(pc_48h,8);	de=ic_48h/100f;
		System.out.println(String.format("48h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_48h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC48M)));
		tmp=countIC(pc_54h);	tmp2=countICAmount(pc_54h,9);	de=ic_54h/100f;
		System.out.println(String.format("54h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_54h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC54M)));
		tmp=countIC(pc_60h);	tmp2=countICAmount(pc_60h,10);	de=ic_60h/100f;
		System.out.println(String.format("60h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_60h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC60M)));
		tmp=countIC(pc_66h);	tmp2=countICAmount(pc_66h,11);	de=ic_66h/100f;
		System.out.println(String.format("66h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_66h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC66M)));
		tmp=countIC(pc_72h);	tmp2=countICAmount(pc_72h,12);	de=ic_72h/100f;
		System.out.println(String.format("72h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_72h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC72M)));
		
		countType(res1,"EFC interaction");
		
		
		System.out.println("\neffective TCs  :"+effectiveTC2);
		System.out.println("effective count:"+effectiveCount2);
		
		ic_6h =NHourIntensityChangeCount(res2,1 );
		ic_12h=NHourIntensityChangeCount(res2,2 );
		ic_18h=NHourIntensityChangeCount(res2,3 );
		ic_24h=NHourIntensityChangeCount(res2,4 );
		ic_30h=NHourIntensityChangeCount(res2,5 );
		ic_36h=NHourIntensityChangeCount(res2,6 );
		ic_42h=NHourIntensityChangeCount(res2,7 );
		ic_48h=NHourIntensityChangeCount(res2,8 );
		ic_54h=NHourIntensityChangeCount(res2,9 );
		ic_60h=NHourIntensityChangeCount(res2,10);
		ic_66h=NHourIntensityChangeCount(res2,11);
		ic_72h=NHourIntensityChangeCount(res2,12);
		
		System.out.println("\ncounts: "+
			ic_6h +"  "+ic_12h+"  "+ic_18h+"  "+ic_24h+"  "+
			ic_30h+"  "+ic_36h+"  "+ic_42h+"  "+ic_48h+"  "+
			ic_54h+"  "+ic_60h+"  "+ic_66h+"  "+ic_72h+"\n"
		);
		
		EFC6M =new float[ic_6h ]; pc_6h =new float[ic_6h ];
		EFC12M=new float[ic_12h]; pc_12h=new float[ic_12h];
		EFC18M=new float[ic_18h]; pc_18h=new float[ic_18h];
		EFC24M=new float[ic_24h]; pc_24h=new float[ic_24h];
		EFC30M=new float[ic_30h]; pc_30h=new float[ic_30h];
		EFC36M=new float[ic_36h]; pc_36h=new float[ic_36h];
		EFC42M=new float[ic_42h]; pc_42h=new float[ic_42h];
		EFC48M=new float[ic_48h]; pc_48h=new float[ic_48h];
		EFC54M=new float[ic_54h]; pc_54h=new float[ic_54h];
		EFC60M=new float[ic_60h]; pc_60h=new float[ic_60h];
		EFC66M=new float[ic_66h]; pc_66h=new float[ic_66h];
		EFC72M=new float[ic_72h]; pc_72h=new float[ic_72h];
		
		pt_6h =0;pt_12h=0;pt_18h=0;pt_24h=0;
		pt_30h=0;pt_36h=0;pt_42h=0;pt_48h=0;
		pt_54h=0;pt_60h=0;pt_66h=0;pt_72h=0;
		for(InteractionEvent e:res2){
			if(e.hasIC(1 )){ EFC6M [pt_6h ]=e.cEFCMean(); pc_6h [pt_6h ]=e.IC(1 ); pt_6h ++;}
			if(e.hasIC(2 )){ EFC12M[pt_12h]=e.cEFCMean(); pc_12h[pt_12h]=e.IC(2 ); pt_12h++;}
			if(e.hasIC(3 )){ EFC18M[pt_18h]=e.cEFCMean(); pc_18h[pt_18h]=e.IC(3 ); pt_18h++;}
			if(e.hasIC(4 )){ EFC24M[pt_24h]=e.cEFCMean(); pc_24h[pt_24h]=e.IC(4 ); pt_24h++;}
			if(e.hasIC(5 )){ EFC30M[pt_30h]=e.cEFCMean(); pc_30h[pt_30h]=e.IC(5 ); pt_30h++;}
			if(e.hasIC(6 )){ EFC36M[pt_36h]=e.cEFCMean(); pc_36h[pt_36h]=e.IC(6 ); pt_36h++;}
			if(e.hasIC(7 )){ EFC42M[pt_42h]=e.cEFCMean(); pc_42h[pt_42h]=e.IC(7 ); pt_42h++;}
			if(e.hasIC(8 )){ EFC48M[pt_48h]=e.cEFCMean(); pc_48h[pt_48h]=e.IC(8 ); pt_48h++;}
			if(e.hasIC(9 )){ EFC54M[pt_54h]=e.cEFCMean(); pc_54h[pt_54h]=e.IC(9 ); pt_54h++;}
			if(e.hasIC(10)){ EFC60M[pt_60h]=e.cEFCMean(); pc_60h[pt_60h]=e.IC(10); pt_60h++;}
			if(e.hasIC(11)){ EFC66M[pt_66h]=e.cEFCMean(); pc_66h[pt_66h]=e.IC(11); pt_66h++;}
			if(e.hasIC(12)){ EFC72M[pt_72h]=e.cEFCMean(); pc_72h[pt_72h]=e.IC(12); pt_72h++;}
		}
		
		de=0;
		tmp=countIC(pc_6h);		tmp2=countICAmount(pc_6h,1);		de=ic_6h /100f;
		System.out.println(String.format(" 6h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_6h ,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC6M)));
		tmp=countIC(pc_12h);	tmp2=countICAmount(pc_12h,2);	de=ic_12h/100f;
		System.out.println(String.format("12h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_12h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC12M)));
		tmp=countIC(pc_18h);	tmp2=countICAmount(pc_18h,3);	de=ic_18h/100f;
		System.out.println(String.format("18h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_18h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC18M)));
		tmp=countIC(pc_24h);	tmp2=countICAmount(pc_24h,4);	de=ic_24h/100f;
		System.out.println(String.format("24h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_24h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC24M)));
		tmp=countIC(pc_30h);	tmp2=countICAmount(pc_30h,5);	de=ic_30h/100f;
		System.out.println(String.format("30h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_30h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC36M)));
		tmp=countIC(pc_36h);	tmp2=countICAmount(pc_36h,6);	de=ic_36h/100f;
		System.out.println(String.format("36h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_36h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC48M)));
		tmp=countIC(pc_42h);	tmp2=countICAmount(pc_42h,7);	de=ic_42h/100f;
		System.out.println(String.format("42h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_42h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC54M)));
		tmp=countIC(pc_48h);	tmp2=countICAmount(pc_48h,8);	de=ic_48h/100f;
		System.out.println(String.format("48h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_48h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC60M)));
		tmp=countIC(pc_54h);	tmp2=countICAmount(pc_54h,9);	de=ic_54h/100f;
		System.out.println(String.format("54h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_54h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC66M)));
		tmp=countIC(pc_60h);	tmp2=countICAmount(pc_60h,10);	de=ic_60h/100f;
		System.out.println(String.format("60h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_60h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC72M)));
		tmp=countIC(pc_66h);	tmp2=countICAmount(pc_66h,11);	de=ic_66h/100f;
		System.out.println(String.format("66h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_66h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC6M)));
		tmp=countIC(pc_72h);	tmp2=countICAmount(pc_72h,12);	de=ic_72h/100f;
		System.out.println(String.format("72h (%4d): %4.1f (%4.1f,%4.1f,%4.1f), %4.1f, %4.1f (MEFC: %5.2f)",pt_72h,
			tmp[0]/de,tmp2[0]/de,tmp2[1]/de,tmp2[2]/de,tmp[1]/de,tmp[2]/de,StatisticsUtil.cArithmeticMean(EFC6M)));
		
		countType(res2,"FFC interaction");
	}
	
	static int[] countIC(float[] PC){
		int[] re=new int[3];
		
		for(int l=0,L=PC.length;l<L;l++)
		if(PC[l]<0) re[0]++;
		else if(PC[l]==0) re[1]++;
		else re[2]++;
		
		return re;
	}
	
	static int[] countICAmount(float[] PC,int h){
		int[] re=new int[3];
		
		for(int l=0,L=PC.length;l<L;l++){
			float tmp=PC[l]/h;
			
			if(tmp<0&&tmp>-5) re[0]++;
			else if(tmp<=-5&&tmp>-10) re[1]++;
			else if(tmp<=-10) re[2]++;
		}
		
		return re;
	}
	
	static InteractionEvent[] getEvents(Typhoon tr,boolean[] validEFC,boolean[] validPCH,boolean[] landing,
	boolean[] hasPoten,float[] EFC,float[] SST,float[] VWS,float[] POT){
		final int continuousSamples=3;
		
		int[] con=getContinuous(validEFC);
		
		InteractionEvent[] res=null;
		
		int size=0;
		
		for(int i=0,str=0;i<con.length;i++){
			if(con[i]>=continuousSamples&&validEFC[str]) size++;
			str+=con[i];
		}
		
		if(size>0){
			res=new InteractionEvent[size];
			
			for(int i=0,str=0,tag=0;i<con.length;i++){
				if(con[i]>=continuousSamples&&validEFC[str]){
					res[tag]=new InteractionEvent(tr,validEFC,validPCH,landing,hasPoten,str,con[i]);
					res[tag].setData(EFC,SST,VWS,POT);
					tag++;
				}
				
				if(tag==size) break;
				
				str+=con[i];
			}
		}
		
		return res;
	}
	
	static int[] getContinuous(boolean[] valid){
		int count=0;
		
		boolean tmp=valid[0];
		
		for(int l=1;l<valid.length;l++)
		if(valid[l]!=tmp){
			count++;
			tmp=valid[l];
		}
		
		int[] con=new int[count+1];
		tmp=valid[0];
		con[0]++;
		
		for(int l=1,tag1=0;l<valid.length;l++){
			if(valid[l]!=tmp){
				tmp=valid[l];
				tag1++;
			}
			
			con[tag1]++;
		}
		
		return con;
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
		System.out.println("excluding littlePoten: "+hasPotent);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
	}
	
	static int NHourIntensityChangeCount(List<InteractionEvent> ls,int offset){
		int count=0;
		
		for(InteractionEvent e:ls)
		if(e.hasIC(offset)) count++;
		
		return count;
	}
	
	static void countType(List<InteractionEvent> ls,String info){
		int[] count=new int[5]; // TD, TS, TY, EC, Others
		
		for(InteractionEvent e:ls){
			int str=e.str;
			TYPE t =e.tr.getTypes()[str];
			
			switch(t){
			case TD:
				count[0]++;	// TD
				break;
			case TS:
				count[1]++;	// TS
				break;
			case TY:
				count[2]++;	// TY
				break;
			case EC:
				count[3]++;	// EC
				break;
			case OTHERS:
				count[4]++;	// OTHERS
				break;

			default:
				throw new IllegalArgumentException("unknown type");
			}
			
			//if(t==TYPE.EC) System.out.println(e);
		}
		
		System.out.println(String.format(
			"\n"+info+": TD(%3d), TS(%3d), TY(%3d), EC(%3d)",
			count[0],count[1],count[2],count[3]
		));
	}
	
	
	@SuppressWarnings("unchecked")
	static List<Event>[] groupEvents(Event[] events){
		List<Event> ls1=new ArrayList<>();
		List<Event> ls2=new ArrayList<>();
		
		for(Event e:events){
			if(e.pslope<=0) ls1.add(e);
			else ls2.add(e);
		}
		
		return new List[]{ls1,ls2};
	}
	
	static void compositeEvents(List<Event> events){
		int len=events.size();
		
		float mEFC=0;
		float mSST=0;
		float mVWS=0;
		float mPRS=0;
		float mSLP=0;
		
		for(Event e:events){
			mEFC+=e.mefc;
			mSST+=e.msst;
			mVWS+=e.mvws;
			mPRS+=e.mprs;
			mSLP+=e.pslope;
		}
		
		mEFC/=len;
		mSST/=len;
		mVWS/=len;
		mPRS/=len;
		mSLP/=len;
		
		System.out.println(String.format("Samples(%4d) EFC(%7.2f) SST(%7.2f) VWS(%7.2f) PRS(%7.2f) SLOPE(%7.2f)",len,mEFC,mSST,mVWS,mPRS,mSLP));
	}
	
	
	/**
	 * class for interaction event
	 */
	static final class InteractionEvent{
		//
		private int N  =0;
		private int str=0;
		
		private float[] EFC=null;
		private float[] SST=null;
		private float[] VWS=null;
		private float[] POT=null;
		private float[] pch=null;
		
		private boolean[] hasPoten=null;
		private boolean[] landing =null;
		private boolean[] validEFC=null;
		private boolean[] validPCH=null;
		
		private Typhoon tr=null;
		
		
		// contructor
		public InteractionEvent(Typhoon tr,boolean[] validEFC,boolean[] validPCH,boolean[] landing,
			boolean[] hasPoten,int l,int len){
			N=len;	this.validEFC=validEFC;	this.landing =landing;
			str=l;	this.validPCH=validPCH;	this.hasPoten=hasPoten;	this.tr=tr;
			
			
			EFC=new float[N];	SST=new float[N];
			VWS=new float[N];	pch=new float[N];	POT=new float[N];
		}
		
		// getor and setor
		public void setData(float[] EFC,float[] SST,float[] VWS,float[] POT){
			System.arraycopy(EFC,str,this.EFC,0,N);
			System.arraycopy(SST,str,this.SST,0,N);
			System.arraycopy(VWS,str,this.VWS,0,N);
			System.arraycopy(POT,str,this.POT,0,N);
			
			float[] pres=tr.getPressures();
			for(int l=0;l<N;l++)
			pch[l]=pres[str+l]-pres[str];
		}
		
		public boolean[] getValidEFC(){ return validEFC;}
		
		public boolean[] getValidPCH(){ return validPCH;}
		
		public boolean[] getLanding(){ return landing;}
		
		public boolean[] getHasPotential(){ return hasPoten;}
		
		public boolean hasIC(int offset){
			if(offset<=0) throw new IllegalArgumentException("offset should be larger than 0");
			
			if(offset+str>=tr.getTCount()) return false;
			
			if(!validPCH[offset+str]) return false;
			
			return true;
		}
		
		public float IC(int offset){
			float[] pres=tr.getPressures();
			
			return pres[offset+str]-pres[str];
		}
		
		public float cEFCMean(){ return StatisticsUtil.cArithmeticMean(EFC);}
		
		
		// print results
		public String toString(){
			StringBuilder sb=new StringBuilder();
			
			sb.append(String.format("%8s(%4s)   Lats     Lons     EFC   FFC   SST   VWS    POT       PRS    PCH  landing  hasPOT\n",
				tr.getName(),tr.getID()
			));
			
			long[] times=tr.getTimes();
			float[] pres=tr.getPressures();
			
			for(int l=str,L=str+N;l<L;l++)
			sb.append(String.format("%10d, %4.1f��N, %5.1f��E, %5.1f, %5.1f, %4.1f, %4.1f, %5.1f, %7.1f, %5.1f,   %5b, %5b\n",
				times[l],tr.getYPositions()[l],tr.getXPositions()[l],
				EFC[l-str],SST[l-str],VWS[l-str],POT[l-str],pres[l],pch[l-str],landing[l],hasPoten[l]
			));
			
			return sb.toString();
		}
	}
	
	static final class Event{
		//
		private int len=0;
		
		private float pslope=0;
		private float eslope=0;
		private float sslope=0;
		private float vslope=0;
		
		private float mprs=0;
		private float mefc=0;
		private float msst=0;
		private float mvws=0;
		
		private float[] EFC=null;
		private float[] SST=null;
		private float[] VWS=null;
		private float[] prs=null;
		
		public Event(InteractionEvent e){
			this.len=e.N;
			this.EFC=new float[len];
			this.SST=new float[len];
			this.VWS=new float[len];
			this.prs=new float[len];
			
			float[] pres=e.tr.getPressures();
			
			for(int l=0,str=e.str;l<len;l++){
				EFC[l]=e.EFC[l];
				SST[l]=e.SST[l];
				VWS[l]=e.VWS[l];
				prs[l]=pres[l+str];
			}
			
			pslope=FilterModel.getLinearTrend(prs);
			eslope=FilterModel.getLinearTrend(EFC);
			sslope=FilterModel.getLinearTrend(SST);
			vslope=FilterModel.getLinearTrend(VWS);
			
			mprs=StatisticsUtil.cArithmeticMean(prs);
			mefc=StatisticsUtil.cArithmeticMean(EFC);
			msst=StatisticsUtil.cArithmeticMean(SST);
			mvws=StatisticsUtil.cArithmeticMean(VWS);
		}
		
		public float getPSlope(){ return pslope;}
		public float getESlope(){ return eslope;}
		public float getSSlope(){ return sslope;}
		public float getVSlope(){ return vslope;}
		
		public float getMEFC(){ return mefc;}
		public float getMSST(){ return msst;}
		public float getMVWS(){ return mvws;}
		public float getMPRS(){ return mprs;}
		
		public void expandTo(int nlen){
			this.len=nlen;
			
			EFC=InterpolationModel.interp1D(EFC,nlen,Type.LINEAR);
			SST=InterpolationModel.interp1D(SST,nlen,Type.LINEAR);
			VWS=InterpolationModel.interp1D(VWS,nlen,Type.LINEAR);
			prs=InterpolationModel.interp1D(prs,nlen,Type.LINEAR);
			
			pslope=FilterModel.getLinearTrend(prs);
			eslope=FilterModel.getLinearTrend(EFC);
			sslope=FilterModel.getLinearTrend(SST);
			vslope=FilterModel.getLinearTrend(VWS);
			
			mprs=StatisticsUtil.cArithmeticMean(prs);
			mefc=StatisticsUtil.cArithmeticMean(EFC);
			msst=StatisticsUtil.cArithmeticMean(SST);
			mvws=StatisticsUtil.cArithmeticMean(VWS);
		}
	}
}
