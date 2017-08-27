package Package;
//
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
import miniufo.statistics.StatisticsUtil;
import static miniufo.statistics.StatisticsUtil.cCorrelationCoefficient;
import static miniufo.test.statistics.SignificanceTest.testCorrelationCoefficient;

//
public class Correlation_bak{
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
	//private static final String respath="d:/Data/PhD/Climatology/"+dsets+"/";
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
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
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
				IntensityModel.greaterEqualThan(tr.getXPositions(),100):IntensityModel.newBooleans(tr.getTCount());
			boolean[] rlon=withinWNP?	// lons <= 190E
				IntensityModel.lessEqualThan(tr.getXPositions(),190):IntensityModel.newBooleans(tr.getTCount());
			
			float[] EFC=efcsm.getData()[1][0][0];
			float[] VWS= vwsm.getData()[0][0][0];
			float[] SST= sstm.getData()[0][0][0];
			float[] POT=IntensityModel.cPotential(SST,tr.getWinds());
			
			boolean[] hEFC=IntensityModel.greaterEqualThan(EFC,10);
			boolean[] hPOT=IntensityModel.greaterEqualThan(POT,30);
			
			boolean[] validEFC=IntensityModel.combination(wind,lsmb,llon,rlon,hEFC,hPOT);
			boolean[] validPCH=IntensityModel.combination(wind,lsmb,llon,rlon,hPOT);
			
			/**
			for(int l=0;l<tr.getCount();l++)
			System.out.println(String.format("%10d, %4.1f¡ãN, %5.1f¡ãE, %5.0f, %5.1f, %4.1f, %4.1f, %4.1f, %6s, %6s",
				tr.getDates()[l].getLongTime(),tr.getLats()[l],tr.getLons()[l],tr.getPressures()[l],
				EFC[l],VWS[l],SST[l],POT[l],validEFC[l],validPCH[l]
			));*/
			
			Event[] events=getEvents(tr,validEFC,validPCH,EFC,SST,VWS,POT);
			
			if(events!=null)
			for(Event e:events){ res.add(e); System.out.println(e);}
			
			int cc=0;
			if(events!=null) cc=events.length;
			effectiveCount+=cc;
			if(cc>0) effectiveTC++;
		}
		
		printResult();
		
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
		
		float[] EFC6M =new float[ic_6h ]; float[] pc_6h =new float[ic_6h ];
		float[] EFC12M=new float[ic_12h]; float[] pc_12h=new float[ic_12h];
		float[] EFC18M=new float[ic_18h]; float[] pc_18h=new float[ic_18h];
		float[] EFC24M=new float[ic_24h]; float[] pc_24h=new float[ic_24h];
		float[] EFC30M=new float[ic_30h]; float[] pc_30h=new float[ic_30h];
		float[] EFC36M=new float[ic_36h]; float[] pc_36h=new float[ic_36h];
		float[] EFC42M=new float[ic_42h]; float[] pc_42h=new float[ic_42h];
		float[] EFC48M=new float[ic_48h]; float[] pc_48h=new float[ic_48h];
		
		int pt_6h =0;
		int pt_12h=0;
		int pt_18h=0;
		int pt_24h=0;
		int pt_30h=0;
		int pt_36h=0;
		int pt_42h=0;
		int pt_48h=0;
		for(Event e:res){
			if(e.hasIntensityChange(1)){ EFC6M [pt_6h ]=e.cEFCMean(); pc_6h [pt_6h ]=e.cIntensityChange(1); pt_6h ++;}
			if(e.hasIntensityChange(2)){ EFC12M[pt_12h]=e.cEFCMean(); pc_12h[pt_12h]=e.cIntensityChange(2); pt_12h++;}
			if(e.hasIntensityChange(3)){ EFC18M[pt_18h]=e.cEFCMean(); pc_18h[pt_18h]=e.cIntensityChange(3); pt_18h++;}
			if(e.hasIntensityChange(4)){ EFC24M[pt_24h]=e.cEFCMean(); pc_24h[pt_24h]=e.cIntensityChange(4); pt_24h++;}
			if(e.hasIntensityChange(5)){ EFC30M[pt_30h]=e.cEFCMean(); pc_30h[pt_30h]=e.cIntensityChange(5); pt_30h++;}
			if(e.hasIntensityChange(6)){ EFC36M[pt_36h]=e.cEFCMean(); pc_36h[pt_36h]=e.cIntensityChange(6); pt_36h++;}
			if(e.hasIntensityChange(7)){ EFC42M[pt_42h]=e.cEFCMean(); pc_42h[pt_42h]=e.cIntensityChange(7); pt_42h++;}
			if(e.hasIntensityChange(8)){ EFC48M[pt_48h]=e.cEFCMean(); pc_48h[pt_48h]=e.cIntensityChange(8); pt_48h++;}			
		}
		
		System.out.println("EFC-6h  IC corr: "+cCorrelationCoefficient(EFC6M ,pc_6h )+" ("+pt_6h +") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC6M ,pc_6h ),0.95f,pt_6h ));
		System.out.println("EFC-12h IC corr: "+cCorrelationCoefficient(EFC12M,pc_12h)+" ("+pt_12h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC12M,pc_12h),0.95f,pt_12h));
		System.out.println("EFC-18h IC corr: "+cCorrelationCoefficient(EFC18M,pc_18h)+" ("+pt_18h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC18M,pc_18h),0.95f,pt_18h));
		System.out.println("EFC-24h IC corr: "+cCorrelationCoefficient(EFC24M,pc_24h)+" ("+pt_24h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC24M,pc_24h),0.95f,pt_24h));
		System.out.println("EFC-30h IC corr: "+cCorrelationCoefficient(EFC30M,pc_30h)+" ("+pt_30h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC30M,pc_30h),0.95f,pt_30h));
		System.out.println("EFC-36h IC corr: "+cCorrelationCoefficient(EFC36M,pc_36h)+" ("+pt_36h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC36M,pc_36h),0.95f,pt_36h));
		System.out.println("EFC-42h IC corr: "+cCorrelationCoefficient(EFC42M,pc_42h)+" ("+pt_42h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC42M,pc_42h),0.95f,pt_42h));
		System.out.println("EFC-48h IC corr: "+cCorrelationCoefficient(EFC48M,pc_48h)+" ("+pt_48h+") "+testCorrelationCoefficient(cCorrelationCoefficient(EFC48M,pc_48h),0.95f,pt_48h));
	}
	
	static Event[] getEvents(Typhoon tr,boolean[] validEFC,boolean[] validPCH,
	float[] EFC,float[] SST,float[] VWS,float[] POT){
		final int continuousSamples=2;
		
		int[] con=getContinuous(validEFC);
		
		Event[] res=null;
		
		int size=0;
		
		for(int i=0,str=0;i<con.length;i++){
			if(con[i]>=continuousSamples&&validEFC[str]) size++;
			str+=con[i];
		}
		
		if(size>0){
			res=new Event[size];
			
			for(int i=0,str=0,tag=0;i<con.length;i++){
				if(con[i]>=continuousSamples&&validEFC[str]){
					res[tag]=new Event(tr,validEFC,validPCH,str,con[i]);
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
	
	static void printResult(){
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
	
	static int NHourIntensityChangeCount(List<Event> ls,int offset){
		int count=0;
		
		for(Event e:ls)
		if(e.hasIntensityChange(offset)) count++;
		
		return count;
	}
	
	/**
	 * class for interaction event
	 */
	static final class Event{
		//
		private int N  =0;
		private int str=0;
		
		private float[] EFC=null;
		private float[] SST=null;
		private float[] VWS=null;
		private float[] POT=null;
		private float[] pch=null;
		
		private boolean[] validEFC=null;
		private boolean[] validPCH=null;
		
		private Typhoon tr=null;
		
		
		// contructor
		public Event(Typhoon tr,boolean[] validEFC,boolean[] validPCH,int l,int len){
			N=len;	this.validEFC=validEFC;
			str=l;	this.validPCH=validPCH;	this.tr=tr;
			
			
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
		
		public boolean hasIntensityChange(int offset){
			if(offset<=0) throw new IllegalArgumentException("offset should be larger than 0");
			
			if(offset+str>=tr.getTCount()) return false;
			
			if(!validPCH[offset+str]) return false;
			
			return true;
		}
		
		public float cIntensityChange(int offset){
			float[] pres=tr.getPressures();
			
			return pres[offset+str]-pres[str];
		}
		
		public float cEFCMean(){ return StatisticsUtil.cArithmeticMean(EFC);}
		
		
		// print results
		public String toString(){
			StringBuilder sb=new StringBuilder();
			
			sb.append(String.format("%8s(%4s)  Lats    Lons      EFC   SST   VWS   POT   PCH\n",
				tr.getName(),tr.getID()
			));
			
			long[] times=tr.getTimes();
			
			for(int l=str,L=str+N;l<L;l++)
			sb.append(String.format("%10d, %4.1f¡ãN, %5.1f¡ãE, %5.1f, %4.1f, %4.1f, %4.1f, %4.1f\n",
				times[l],tr.getYPositions()[l],tr.getXPositions()[l],
				EFC[l-str],SST[l-str],VWS[l-str],POT[l-str],pch[l-str]
			));
			
			return sb.toString();
		}
	}
}
