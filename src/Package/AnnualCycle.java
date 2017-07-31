package Package;
//
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
import miniufo.lagrangian.Typhoon;

//
public class AnnualCycle{
	//
	private static int ystr=1989;
	private static int yend=2009;
	
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
		
		Statis stsALL =new Statis("all samples" );
		Statis stsEFCS=new Statis("|EFCS| >= 10");
		Statis stsEFCL=new Statis("|EFCL| >= 10");
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,25,2,0.3f,-650));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v");
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getZonalVelocity(),tr.getMeridionalVelocity(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efclm=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),15,24);	// 500-800 km
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
			
			boolean[] EFCS=IntensityModel.ABSgreaterEqualThan(efcsm.getData()[1][0][0],10);// |EFCS| >= 10
			boolean[] EFCL=IntensityModel.ABSgreaterEqualThan(efclm.getData()[1][0][0],10);// |EFCL| >= 10
			
			boolean[] validAll =IntensityModel.combination(wind,lsmb,llon,rlon);
			boolean[] validEFCS=IntensityModel.combination(validAll,EFCS);
			boolean[] validEFCL=IntensityModel.combination(validAll,EFCL);
			
			int cc=IntensityModel.getValidCount(validAll);
			effectiveCount+=cc;
			if(cc>0) effectiveTC++;
			if(effectiveTC%20==0) System.out.print(".");
			
			 stsALL.compile(tr,validAll );
			stsEFCS.compile(tr,validEFCS);
			stsEFCL.compile(tr,validEFCL);
		}
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
		
		System.out.println( stsALL);
		System.out.println(stsEFCS);
		System.out.println(stsEFCL);
	}
	
	
	//
	static final class Statis{
		//
		private int[] splCC=null;
		private int[] TC_CC=null;
		
		private String name =null;
		
		
		//
		public Statis(String name){
			this.name=name;
			
			splCC=new int[12];
			TC_CC=new int[12];
		}
		
		public void compile(Typhoon tr,boolean[] valid){
			int len=tr.getTCount();
			int vc =IntensityModel.getValidCount(valid);
			
			if(len!=valid.length)
			throw new IllegalArgumentException("lengths not equal");
			
			long[] times=tr.getTimes();
			float[] wnds=tr.getWinds();
			
			for(int l=0;l<len;l++)
			if(vc>0&&wnds[l]>17.2f){
				switch(getMonth(times[l])){
					case  1:	TC_CC[ 0]++; break;
					case  2:	TC_CC[ 1]++; break;
					case  3:	TC_CC[ 2]++; break;
					case  4:	TC_CC[ 3]++; break;
					case  5:	TC_CC[ 4]++; break;
					case  6:	TC_CC[ 5]++; break;
					case  7:	TC_CC[ 6]++; break;
					case  8:	TC_CC[ 7]++; break;
					case  9:	TC_CC[ 8]++; break;
					case 10:	TC_CC[ 9]++; break;
					case 11:	TC_CC[10]++; break;
					case 12:	TC_CC[11]++; break;
					default: throw new IllegalArgumentException("invalid month");
				}
				
				break;
			}
			
			for(int l=0;l<len;l++)
			if(valid[l])
			switch(getMonth(times[l])){
				case  1:	splCC[ 0]++; break;
				case  2:	splCC[ 1]++; break;
				case  3:	splCC[ 2]++; break;
				case  4:	splCC[ 3]++; break;
				case  5:	splCC[ 4]++; break;
				case  6:	splCC[ 5]++; break;
				case  7:	splCC[ 6]++; break;
				case  8:	splCC[ 7]++; break;
				case  9:	splCC[ 8]++; break;
				case 10:	splCC[ 9]++; break;
				case 11:	splCC[10]++; break;
				case 12:	splCC[11]++; break;
				default: throw new IllegalArgumentException("invalid month");
			}
		}
		
		public String toString(){
			int splAll=0;	for(int i:splCC) splAll+=i;
			int TC_All=0;	for(int i:TC_CC) TC_All+=i;
			
			StringBuilder sb=new StringBuilder();
			
			sb.append("\n----------------------------- Annual cycle statistics for "+name+" --------------------------\n");
			sb.append("   month         :   Jan   Feb   Mar   Apr   May   Jun   Jul   Aug   Sep   Oct   Nov   Dec,   All\n");
			sb.append("TC count "+(yend-ystr+1)+" years:");
			for(int i=0;i<12;i++) sb.append(String.format("%6s",TC_CC[i]));
			
			sb.append(","+String.format("%6s",TC_All)+"\n");
			
			sb.append("TC count per year:");
			for(int i=0;i<12;i++) sb.append(String.format("%6.1f",TC_CC[i]/(float)(yend-ystr+1)));
			
			sb.append(","+String.format("%6.1f",TC_All/(float)(yend-ystr+1))+"\n");
			
			sb.append("SP count "+(yend-ystr+1)+" years:");
			for(int i=0;i<12;i++) sb.append(String.format("%6s",splCC[i]));
			
			sb.append(","+String.format("%6s",splAll)+"\n");
			
			sb.append("SP count per year:");
			for(int i=0;i<12;i++) sb.append(String.format("%6.1f",splCC[i]/(float)(yend-ystr+1)));
			
			sb.append(","+String.format("%6.1f",splAll/(float)(yend-ystr+1))+"\n");
			
			sb.append("-------------------------------------------------------------------------------------------------\n");
			
			return sb.toString();
		}
		
		private static int getMonth(long time){
			return (int)(time%10000000000L/100000000L);
		}
	}
}
