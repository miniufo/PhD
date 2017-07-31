package Package;

import java.io.FileWriter;
import java.io.IOException;
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
import miniufo.diagnosis.Variable.Dimension;
import miniufo.lagrangian.Typhoon;

//
public class StatFactors{
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
		
		Stat vwsstat =new Stat(new float[]{0,5,10,15,                    Float.MAX_VALUE},"VWS" ,respath);
		Stat sststat =new Stat(new float[]{25.5f,26.5f,27.5f,28.5f,29.5f,Float.MAX_VALUE},"SST" ,respath);
		Stat mpistat =new Stat(new float[]{0,10,20,30,40,50,60,70,80,90, Float.MAX_VALUE},"MPI" ,respath);
		Stat efclstat=new Stat(new float[]{-20,-10,0,10,20,              Float.MAX_VALUE},"EFCL",respath);
		Stat efcsstat=new Stat(new float[]{-20,-10,0,10,20,              Float.MAX_VALUE},"EFCS",respath);
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,28,2,0.3f,-650,850));
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
			Variable efclm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y,15,24);	// 500-800 km
			Variable efcsm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y, 9,18);	// 300-600 km
			
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
			
			boolean[] vali=IntensityModel.combination(wind,lsmb,llon,rlon);
			
			int cc=IntensityModel.getValidCount(vali);
			effectiveCount+=cc;
			if(cc>0){
				effectiveTC++;
				if(effectiveTC%20==0) System.out.print(".");
			}
			
			float[] potential=IntensityModel.cPotential(sstm.getData()[0][0][0],tr.getWinds());
			
			 vwsstat.compileByPressure(tr, vwsm.getData()[0][0][0],vali);
			 sststat.compileByPressure(tr, sstm.getData()[0][0][0],vali);
			 mpistat.compileByPressure(tr,potential,vali);
			efclstat.compileByPressure(tr,efclm.getData()[1][0][0],vali);
			efcsstat.compileByPressure(tr,efcsm.getData()[1][0][0],vali);
		}
		
		printResult(vwsstat,sststat,mpistat,efclstat,efcsstat);
	}
	
	static void printResult(Stat... s){
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
		System.out.println("forward diff for delP: "+forwardDf);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
		
		for(int i=0,I=s.length;i<I;i++){
			System.out.println(s[i]+"\n");
			s[i].release();
		}
	}
	
	
	/**
	 * class for statistics
	 */
	static final class Stat{
		//
		int N=0;
		
		float[] splitor=null;
		
		int[][] count=null;	// [0] is intensify, [1] is no change and [2] is weaken
		
		String name=null;
		String path=null;
		
		FileWriter[][] fw=null;
		
		
		// contructor
		public Stat(float[] sep,String vname,String fpath){
			N=sep.length;
			
			path   =fpath;
			name   =vname;
			splitor=sep;
			count  =new int[3][N];
			fw     =new FileWriter[3][N];
			
			try{
				for(int i=0;i<N;i++){
					fw[0][i]=new FileWriter(path+"TXT/"+name+i+"S.txt");
					fw[1][i]=new FileWriter(path+"TXT/"+name+i+"M.txt");
					fw[2][i]=new FileWriter(path+"TXT/"+name+i+"W.txt");
				}
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		// add up data
		public void compileByPressure(Typhoon tr,float[] data,boolean[] valid){
			int len=tr.getTCount();
			
			if(len!=data.length||len!=valid.length)
			throw new IllegalArgumentException("lengths not equal");
			
			float[] dpr=forwardDf?
				Typhoon.getChangesByForwardDiff(tr.getPressures(),1):
				Typhoon.getChangesByCentralDiff(tr.getPressures());
			float[] lon=tr.getLongitudes();
			float[] lat=tr.getLatitudes();
			
			try{
				for(int l=0;l<len;l++)
				if(valid[l])
				for(int i=0;i<N;i++)
				if(data[l]-splitor[i]<0)
				if(dpr[l]<0){
					count[0][i]++;
					
					fw[0][i].write(lon[l]+" "+lat[l]+"\n");
					
					break;
					
				}else if(dpr[l]==0){
					count[1][i]++;
					
					fw[1][i].write(lon[l]+" "+lat[l]+"\n");
					
					break;
					
				}else{
					count[2][i]++;
					
					fw[2][i].write(lon[l]+" "+lat[l]+"\n");
					
					break;
				}
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		public void compileByWind(Typhoon tr,float[] data,boolean[] valid){
			int len=tr.getTCount();
			
			if(len!=data.length||len!=valid.length)
			throw new IllegalArgumentException("lengths not equal");
			
			float[] dwd=Typhoon.getChangesByCentralDiff(tr.getWinds());
			float[] lon=tr.getLongitudes();
			float[] lat=tr.getLatitudes();
			
			try{
				for(int l=0;l<len;l++)
				if(valid[l])
				for(int i=0;i<N;i++)
				if(data[l]-splitor[i]<0)
				if(dwd[l]>0){
					count[0][i]++;
					
					fw[0][i].write(lon[l]+" "+lat[l]+"\n");
					
					break;
					
				}else if(dwd[l]==0){
					count[1][i]++;
					
					fw[1][i].write(lon[l]+" "+lat[l]+"\n");
					
					break;
					
				}else{
					count[2][i]++;
					
					fw[2][i].write(lon[l]+" "+lat[l]+"\n");
					
					break;
				}
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		// close all files
		public void release(){
			try{
				for(int i=0;i<N;i++){
					fw[0][i].close();
					fw[1][i].close();
					fw[2][i].close();
				}
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
		
		
		// print results
		public String toString(){
			StringBuilder sb=new StringBuilder();
			
			sb.append("----------------statistics of "+name+"----------------\n");
			sb.append("    ranging       total   deepen   remain    fill\n");
			
			sb.append("(   -¡Þ, "+String.format("%5.1f",splitor[0])+"):"+
				String.format("%8d",count[0][0]+count[1][0]+count[2][0])+
				String.format("%9d",count[0][0])+
				String.format("%9d",count[1][0])+
				String.format("%8d",count[2][0])+"\n"
			);
			
			for(int i=1;i<N-1;i++)
			sb.append(
				"["+
				String.format("%5.1f",splitor[i-1])+
				", "+
				String.format("%5.1f",splitor[i])+
				"):"+
				String.format("%8d",count[0][i]+count[1][i]+count[2][i])+
				String.format("%9d",count[0][i])+
				String.format("%9d",count[1][i])+
				String.format("%8d",count[2][i])+"\n"
			);
			
			sb.append("["+String.format("%5.1f",splitor[N-2])+",    +¡Þ):"+
				String.format("%8d",count[0][N-1]+count[1][N-1]+count[2][N-1])+
				String.format("%9d",count[0][N-1])+
				String.format("%9d",count[1][N-1])+
				String.format("%8d",count[2][N-1])+"\n"
			);
			
			sb.append("-------------------------------------------------\n");
			
			return sb.toString();
		}
	}
}
