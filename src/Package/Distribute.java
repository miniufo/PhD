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
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;


//
public class Distribute{
	//
	private static int effectiveCount=0;
	private static int effectiveTC   =0;
	
	private static final boolean noDepress=false;	// wind > 17.2 m/s
	private static final boolean noLanding=false;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	
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
		
		Distribution disSST =new Distribution(ctl,"SST" ,new float[]{25.5f,26.5f,27.5f,28.5f,29.5f,Float.MAX_VALUE});
		Distribution disVWS =new Distribution(ctl,"VWS" ,new float[]{0,5,10,15,                    Float.MAX_VALUE});
		Distribution disEFCS=new Distribution(ctl,"EFCS",new float[]{-20,-10,0,10,20,              Float.MAX_VALUE});
		Distribution disEFCL=new Distribution(ctl,"EFCL",new float[]{-20,-10,0,10,20,              Float.MAX_VALUE});
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,25,2,0.3f,-650));
			CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
			
			dd.setCtlDescriptor(ctl);	df.setPrinting(false);
			CylindricalSpatialModel csm=new CylindricalSpatialModel(dd);
			DynamicMethodsInCC dm=new DynamicMethodsInCC(csm);
			CoordinateTransformation ct=new CoordinateTransformation(new SphericalSpatialModel(ctl),csm);
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v","sst");
			Variable[] shrs=dm.cVerticalWindShear(vars[0],vars[1]);
			
			Variable vwsm=dm.cRadialAverage(shrs[0].hypotenuse(shrs[1]),1,15).anomalizeX();
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
			
			boolean[] valid=IntensityModel.combination(wind,lsmb,llon,rlon);
			
			 disSST.compile(tr, sstm.getData()[0][0][0],valid);
			 disVWS.compile(tr, vwsm.getData()[0][0][0],valid);
			disEFCS.compile(tr,efcsm.getData()[1][0][0],valid);
			disEFCL.compile(tr,efclm.getData()[1][0][0],valid);
			
			int cc=IntensityModel.getValidCount(valid);
			effectiveCount+=cc;
			if(cc>0){
				effectiveTC++;
				if(effectiveTC%20==0) System.out.print(".");
			}
		}
		
		System.out.println("\n");
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
		
		 disSST.average();
		 disVWS.average();
		disEFCS.average();
		disEFCL.average();
		
		DataWrite dw=DataIOFactory.getDataWrite(ctl,respath+"Distribute.dat");
		dw.writeData(ctl,ArrayUtil.concatAll(Variable.class,
			 disSST.getLonDisVariables(), disSST.getLatDisVariables(), disSST.getHorizontalVariables(),
			 disVWS.getLonDisVariables(), disVWS.getLatDisVariables(), disVWS.getHorizontalVariables(),
			disEFCS.getLonDisVariables(),disEFCS.getLatDisVariables(),disEFCS.getHorizontalVariables(),
			disEFCL.getLonDisVariables(),disEFCL.getLatDisVariables(),disEFCL.getHorizontalVariables()
		));	dw.closeFile();
	}
	
	
	//
	static final class Distribution{
		//
		private int count=0;
		
		private float[] separator=null;
		
		private Variable[] lonDis=null;
		private Variable[] latDis=null;
		private Variable[] horDis=null;
		
		private DataDescriptor dd=null;
		
		
		//
		public Distribution(DataDescriptor dd,String name,float[] sep){
			this.dd=dd;	this.separator=sep;	count=sep.length;
			
			lonDis=new Variable[count];
			latDis=new Variable[count];
			horDis=new Variable[count];
			
			lonDis[0]=new Variable(name+"Lon"+0,true,new Range("t(1,1);z(1,1)",dd));
			lonDis[0].setCommentAndUnit("longitudinal distribution of "+name+" (-¡Þ, "+sep[0]+")");
			lonDis[0].setUndef(-9999);
			
			latDis[0]=new Variable(name+"Lat"+0,true,new Range("t(1,1);z(1,1)",dd));
			latDis[0].setCommentAndUnit("latitudinal distribution of "+name+" (-¡Þ, "+sep[0]+")");
			latDis[0].setUndef(-9999);
			
			horDis[0]=new Variable(name+"Hor"+0,true,new Range("t(1,1);z(1,1)",dd));
			horDis[0].setCommentAndUnit("horizontal distribution of "+name+" (-¡Þ, "+sep[0]+")");
			horDis[0].setUndef(-9999);
			
			for(int i=1;i<count;i++){
				lonDis[i]=new Variable(name+"Lon"+i,true,new Range("t(1,1);z(1,1)",dd));
				lonDis[i].setCommentAndUnit("longitudinal distribution of "+name+" ("+sep[i-1]+", "+sep[i]+")");
				lonDis[i].setUndef(-9999);
				
				latDis[i]=new Variable(name+"Lat"+i,true,new Range("t(1,1);z(1,1)",dd));
				latDis[i].setCommentAndUnit("latitudinal distribution of "+name+" ("+sep[i-1]+", "+sep[i]+")");
				latDis[i].setUndef(-9999);
				
				horDis[i]=new Variable(name+"Hor"+i,true,new Range("t(1,1);z(1,1)",dd));
				horDis[i].setCommentAndUnit("horizontal distribution of "+name+" ("+sep[i-1]+", "+sep[i]+")");
				horDis[i].setUndef(-9999);
			}
		}
		
		
		// getor and setor
		public Variable[] getLonDisVariables(){ return lonDis;}
		
		public Variable[] getLatDisVariables(){ return latDis;}
		
		public Variable[] getHorizontalVariables(){ return horDis;}
		
		
		//
		public void compile(Typhoon tr,float[] data,boolean[] valid){
			float[] lons=tr.getLongitudes();
			float[] lats=tr.getLatitudes();
			
			for(int l=0,L=tr.getTCount();l<L;l++)
			if(valid[l])
			for(int m=0;m<count;m++)
			if(data[l]-separator[m]<0){
				for(int j=0,J=lonDis[0].getYCount();j<J;j++)
				lonDis[m].getData()[0][0][j][dd.getXNum(lons[l])]++;
				
				for(int i=0,I=lonDis[0].getXCount();i<I;i++)
				latDis[m].getData()[0][0][dd.getYNum(lats[l])][i]++;
				
				horDis[m].getData()[0][0][dd.getYNum(lats[l])][dd.getXNum(lons[l])]++;
				
				break;
			}
		}
		
		public void average(){
			float dx=dd.getDXDef()[0];
			float dy=dd.getDYDef()[0];
			
			for(int i=0;i<count;i++){
				lonDis[i].divideEq(dx);
				latDis[i].divideEq(dy);
			}
		}
		
		
		public static boolean[] isType(TYPE[] types,TYPE type){
			int N=types.length;
			
			boolean[] b=new boolean[N];
			
			for(int l=0;l<N;l++) b[l]=(types[l]==type);
			
			return b;
		}
	}
}
