package Package;
//
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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
import static miniufo.statistics.StatisticsUtil.cCorrelationCoefficient;
import static miniufo.statistics.StatisticsUtil.cPartialCorrelationCoefficient;
import static miniufo.test.statistics.SignificanceTest.testCorrelationCoefficient;

//
public class Correlation{
	//
	private static int effectiveCount=0;
	private static int effectiveTC   =0;
	
	private static float minLat= 90;
	private static float minLon=360;
	private static float maxLat=-90;
	private static float maxLon=  0;
	
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
		
		List<Event> res=new ArrayList<Event>();
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,25,2,0.3f,-650));
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
			dm.cStormRelativeAziRadVelocity(tr.getZonalVelocity(),tr.getMeridionalVelocity(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			//Variable efclm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),15,24);	// 500-800 km
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
			
			float[] EFC=efcsm.getData()[1][0][0];
			float[] VWS= vwsm.getData()[0][0][0];
			float[] SST= sstm.getData()[0][0][0];
			float[] POT=IntensityModel.cPotential(SST,tr.getWinds());
			
			boolean[] hPOT=IntensityModel.greaterEqualThan(POT,30);
			
			boolean[] validPCH=IntensityModel.combination(wind,lsmb,llon,rlon,hPOT);
			
			/**
			for(int l=0;l<tr.getCount();l++)
			System.out.println(String.format("%10d, %4.1f°„N, %5.1f°„E, %5.0f, %5.1f, %4.1f, %4.1f, %4.1f, %6s",
				tr.getDates()[l].getLongTime(),tr.getLats()[l],tr.getLons()[l],tr.getPressures()[l],
				EFC[l],VWS[l],SST[l],POT[l],validPCH[l]
			));*/
			
			Event[] events=getEvents(tr,validPCH,EFC,SST,VWS,POT);
			
			if(events!=null)
			for(Event e:events){ res.add(e);}
			
			int cc=0;
			if(events!=null) cc=events.length;
			effectiveCount+=cc;
			if(cc>0) effectiveTC++;
		}
		
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
		
		int ic_6h =NHourIntensityChangeCount(res,1);
		int ic_12h=NHourIntensityChangeCount(res,2);
		int ic_18h=NHourIntensityChangeCount(res,3);
		int ic_24h=NHourIntensityChangeCount(res,4);
		int ic_30h=NHourIntensityChangeCount(res,5);
		int ic_36h=NHourIntensityChangeCount(res,6);
		int ic_42h=NHourIntensityChangeCount(res,7);
		int ic_48h=NHourIntensityChangeCount(res,8);
		
		System.out.println("counts: "+
			ic_6h +"  "+ic_12h+"  "+ic_18h+"  "+ic_24h+"  "+
			ic_30h+"  "+ic_36h+"  "+ic_42h+"  "+ic_48h
		);
		
		float[] EFC6M =new float[ic_6h ]; float[] VWS6M =new float[ic_6h ];
		float[] EFC12M=new float[ic_12h]; float[] VWS12M=new float[ic_12h];
		float[] EFC18M=new float[ic_18h]; float[] VWS18M=new float[ic_18h];
		float[] EFC24M=new float[ic_24h]; float[] VWS24M=new float[ic_24h];
		float[] EFC30M=new float[ic_30h]; float[] VWS30M=new float[ic_30h];
		float[] EFC36M=new float[ic_36h]; float[] VWS36M=new float[ic_36h];
		float[] EFC42M=new float[ic_42h]; float[] VWS42M=new float[ic_42h];
		float[] EFC48M=new float[ic_48h]; float[] VWS48M=new float[ic_48h];
		
		float[] POT6M =new float[ic_6h ]; float[] pc_6h =new float[ic_6h ];
		float[] POT12M=new float[ic_12h]; float[] pc_12h=new float[ic_12h];
		float[] POT18M=new float[ic_18h]; float[] pc_18h=new float[ic_18h];
		float[] POT24M=new float[ic_24h]; float[] pc_24h=new float[ic_24h];
		float[] POT30M=new float[ic_30h]; float[] pc_30h=new float[ic_30h];
		float[] POT36M=new float[ic_36h]; float[] pc_36h=new float[ic_36h];
		float[] POT42M=new float[ic_42h]; float[] pc_42h=new float[ic_42h];
		float[] POT48M=new float[ic_48h]; float[] pc_48h=new float[ic_48h];
		
		int pt_6h =0,pt_12h=0,pt_18h=0,pt_24h=0,pt_30h=0,pt_36h=0,pt_42h=0,pt_48h=0;
		for(Event e:res){
			if(e.hasIC(1)){ EFC6M [pt_6h ]=e.getEFC(); VWS6M [pt_6h ]=e.getVWS(); POT6M [pt_6h ]=e.getPOT(); pc_6h [pt_6h ]=e.IC(1); pt_6h ++;}
			if(e.hasIC(2)){ EFC12M[pt_12h]=e.getEFC(); VWS12M[pt_12h]=e.getVWS(); POT12M[pt_12h]=e.getPOT(); pc_12h[pt_12h]=e.IC(2); pt_12h++;}
			if(e.hasIC(3)){ EFC18M[pt_18h]=e.getEFC(); VWS18M[pt_18h]=e.getVWS(); POT18M[pt_18h]=e.getPOT(); pc_18h[pt_18h]=e.IC(3); pt_18h++;}
			if(e.hasIC(4)){ EFC24M[pt_24h]=e.getEFC(); VWS24M[pt_24h]=e.getVWS(); POT24M[pt_24h]=e.getPOT(); pc_24h[pt_24h]=e.IC(4); pt_24h++;}
			if(e.hasIC(5)){ EFC30M[pt_30h]=e.getEFC(); VWS30M[pt_30h]=e.getVWS(); POT30M[pt_30h]=e.getPOT(); pc_30h[pt_30h]=e.IC(5); pt_30h++;}
			if(e.hasIC(6)){ EFC36M[pt_36h]=e.getEFC(); VWS36M[pt_36h]=e.getVWS(); POT36M[pt_36h]=e.getPOT(); pc_36h[pt_36h]=e.IC(6); pt_36h++;}
			if(e.hasIC(7)){ EFC42M[pt_42h]=e.getEFC(); VWS42M[pt_42h]=e.getVWS(); POT42M[pt_42h]=e.getPOT(); pc_42h[pt_42h]=e.IC(7); pt_42h++;}
			if(e.hasIC(8)){ EFC48M[pt_48h]=e.getEFC(); VWS48M[pt_48h]=e.getVWS(); POT48M[pt_48h]=e.getPOT(); pc_48h[pt_48h]=e.IC(8); pt_48h++;}			
		}
		
		printResult(EFC6M ,VWS6M ,pc_6h ,"6h ");
		printResult(EFC12M,VWS12M,pc_12h,"12h");
		printResult(EFC18M,VWS18M,pc_18h,"18h");
		printResult(EFC24M,VWS24M,pc_24h,"24h");
		printResult(EFC30M,VWS30M,pc_30h,"30h");
		printResult(EFC36M,VWS36M,pc_36h,"36h");
		printResult(EFC42M,VWS42M,pc_42h,"42h");
		printResult(EFC48M,VWS48M,pc_48h,"48h");
		
		writeData(EFC6M ,VWS6M ,POT6M ,pc_6h ,respath+"/Correlation/Corr_6h.txt" );
		writeData(EFC12M,VWS12M,POT12M,pc_12h,respath+"/Correlation/Corr_12h.txt");
		writeData(EFC18M,VWS18M,POT18M,pc_18h,respath+"/Correlation/Corr_18h.txt");
		writeData(EFC24M,VWS24M,POT24M,pc_24h,respath+"/Correlation/Corr_24h.txt");
		writeData(EFC30M,VWS30M,POT30M,pc_30h,respath+"/Correlation/Corr_30h.txt");
		writeData(EFC36M,VWS36M,POT36M,pc_36h,respath+"/Correlation/Corr_36h.txt");
		writeData(EFC42M,VWS42M,POT42M,pc_42h,respath+"/Correlation/Corr_42h.txt");
		writeData(EFC48M,VWS48M,POT48M,pc_48h,respath+"/Correlation/Corr_48h.txt");
	}
	
	static void writeData(float[] EFC,float[] VWS,float[] POT,float[] PC,String path){
		try{
			FileWriter fw=new FileWriter(new File(path));
			for(int l=0;l<EFC.length;l++)
			fw.write(String.format("%5.1f, %5.1f, %5.1f, %4.1f\n",EFC[l],VWS[l],POT[l],PC[l]));
			fw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static Event[] getEvents(Typhoon tr,boolean[] validPCH,
	float[] EFC,float[] SST,float[] VWS,float[] POT){
		Event[] res=null;
		
		int size=IntensityModel.getValidCount(validPCH);
		
		if(size>0){
			res=new Event[size];
			
			for(int i=0,tag=0;i<tr.getTCount();i++){
				if(validPCH[i]){
					res[tag]=new Event(tr,validPCH,i);
					res[tag].setData(EFC,SST,VWS,POT);
					tag++;
				}
				
				if(tag==size) break;
			}
		}
		
		return res;
	}
	
	static void printResult(float[] EFC,float[] VWS,float[] PC,String hrs){
		System.out.println(String.format(
			"EFC-"+hrs+"  IC corr: %10.8f (%4d) %5s;\tpartial corr: %10.8f",
			cCorrelationCoefficient(EFC,PC),EFC.length,
			testCorrelationCoefficient(cCorrelationCoefficient(EFC,PC),0.95f,EFC.length ),
			cPartialCorrelationCoefficient(EFC,PC,VWS)
		));
	}
	
	static int NHourIntensityChangeCount(List<Event> ls,int offset){
		int count=0;
		
		for(Event e:ls)
		if(e.hasIC(offset)) count++;
		
		return count;
	}
	
	/**
	 * class for interaction event
	 */
	static final class Event{
		//
		private int str=0;
		
		private float EFC;
		private float SST;
		private float VWS;
		private float POT;
		
		private boolean[] validPCH=null;
		
		private Typhoon tr=null;
		
		
		// contructor
		public Event(Typhoon tr,boolean[] validPCH,int l){
			str=l;
			this.tr=tr;
			this.validPCH=validPCH;
		}
		
		// getor and setor
		public void setData(float[] EFC,float[] SST,float[] VWS,float[] POT){
			this.EFC=EFC[str];
			this.SST=SST[str];
			this.VWS=VWS[str];
			this.POT=POT[str];
		}
		
		public boolean[] getValidPCH(){ return validPCH;}
		
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
		
		public float getEFC(){ return EFC;}
		
		public float getSST(){ return SST;}
		
		public float getVWS(){ return VWS;}
		
		public float getPOT(){ return POT;}
		
		
		// print results
		public String toString(){
			StringBuilder sb=new StringBuilder();
			
			sb.append(String.format("%8s(%4s)  Lats    Lons      EFC   SST   VWS   POT   PCH\n",
				tr.getName(),tr.getID()
			));
			
			sb.append(String.format("%10d, %4.1f°„N, %5.1f°„E, %5.1f, %4.1f, %4.1f, %4.1f\n",
				tr.getTime(str),tr.getLatitudes()[str],tr.getLongitudes()[str],EFC,SST,VWS,POT
			));
			
			return sb.toString();
		}
	}
}
