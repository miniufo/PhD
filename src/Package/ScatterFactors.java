package Package;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.basic.ArrayUtil;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.database.DataBaseUtil;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;

//
public class ScatterFactors{
	//
	private static int allCount=0;
	private static int allTC   =0;
	private static int deltaP_interval=1;		// 1 for 6hr, 2 for 12hr...
	
	private static float minLat= 90;
	private static float minLon=360;
	private static float maxLat=-90;
	private static float maxLon=  0;
	
	private static final float undef=DataBaseUtil.undef;
	
	private static final boolean noDepress=true;	// wind > 17.2 m/s
	private static final boolean noLanding=true;	// no land within 200 km
	private static final boolean withinWNP=false;	// 100 <= longitudes <= 190
	private static final boolean forwardDf=true;	// using forward difference to get deltaP
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Statistics/"+dsets+"/Scatter/";
	private static final String tranges="time=1Jan1987-31Dec2011";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		DiagnosisFactory df2=DiagnosisFactory.parseFile("d:/Data/ERAInterim/lsm.ctl");
		
		CtlDescriptor ctl2=(CtlDescriptor)df2.getDataDescriptor();
		Variable lsm=df2.getVariables(new Range("z(1,1)",ctl2),false,"lsm")[0];
		
		for(Typhoon tr:ls){ allTC++; allCount+=tr.getTCount();}
		
		Scatter sctr=new Scatter(allCount,new float[]{Float.NEGATIVE_INFINITY,-50,-30,-20,-10,-5,0,5,10,20,30,50,Float.MAX_VALUE});
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,25,2,0.3f,-650,850));
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
			
			Variable mag =vars[0].hypotenuse(vars[1]);
			Variable magm=dm.cRadialAverage(mag,1,15).anomalizeX();
			
			Variable shrsum=dm.cRadialAverage(shrs[0],1,15).anomalizeX();
			Variable shrsvm=dm.cRadialAverage(shrs[1],1,15).anomalizeX();
			
			Variable vwsm=shrsum.hypotenuse(shrsvm);
			Variable sstm=dm.cRadialAverage(vars[2],1,15).anomalizeX().minusEq(273.15f);
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
			Variable utm=utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable efclm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y,15,24);	// 500-800 km
			Variable efcsm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y, 9,18);	// 300-600 km
			Variable isbym=dm.cMeanInertialStabilityByUT(utm).averageAlong(Dimension.Y,9,18);
			
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
			
			boolean[] vali=IntensityModel.combination(wind,lsmb,llon,rlon);
			
			float[] mwnd=tr.getWinds();
			float[] srwd=tr.getStormRelativeWinds();
			float[] pres=tr.getPressures();
			
			sctr.compile(vali,
				sstm.getData()[0][0][0],												// SST
				vwsm.getData()[0][0][0],												// VWS
				IntensityModel.cMPI(sstm.getData()[0][0][0]),							// MPI
				efcsm.getData()[1][0][0],												// EFCS
				efclm.getData()[1][0][0],												// EFCL
				mwnd,																	// wind
				srwd,																	// storm-relative wind
				pres,																	// pressure
				Typhoon.getChangesByCentralDiff(sstm.getData()[0][0][0]),						// delta SST
				Typhoon.getChangesByCentralDiff(vwsm.getData()[0][0][0]),						// delta VWS
				Typhoon.getChangesByCentralDiff(IntensityModel.cMPI(sstm.getData()[0][0][0])),	// delta MPI
				Typhoon.getChangesByCentralDiff(efcsm.getData()[1][0][0]),						// delta EFCS
				Typhoon.getChangesByCentralDiff(efclm.getData()[1][0][0]),						// delta EFCL
				forwardDf?
				Typhoon.getChangesByForwardDiff(mwnd,deltaP_interval):
				Typhoon.getChangesByCentralDiff(mwnd),									// delta wind
				forwardDf?
				Typhoon.getChangesByForwardDiff(srwd,deltaP_interval):
				Typhoon.getChangesByCentralDiff(srwd),									// delta storm-relative wind
				forwardDf?
				Typhoon.getChangesByForwardDiff(pres,deltaP_interval):					// delta pressure
				Typhoon.getChangesByCentralDiff(pres),						
				tr.getUVel(),													// zonal translating speed
				tr.getVVel(),												// meridional translating speed
				tr.getSpeeds(),															// translating speed
				tr.getXPositions(),														// longitudes
				tr.getYPositions(),														// latitudes
				magm.getData()[1][0][0],
				getTypes(tr),
				isbym.getData()[1][0][0]												// inertial stability
			);
		}
		
		System.out.println("\n\npointer: "+sctr.getPointer());
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("forward diff for delP: "+forwardDf);
		System.out.println("all TCs  :"+allTC);
		System.out.println("all count:"+allCount+"\n");
		
		sctr.count();
		
		sctr.writeToFile(respath+"Scatter"+deltaP_interval+".dat");
	}
	
	private static float[] getTypes(Typhoon tr){
		int count=tr.getTCount();
		
		float[] typeN=new float[count];
		
		TYPE[] types=tr.getTypes();
		
		for(int l=0;l<count;l++){
			switch(types[l]){
			case TD:		typeN[l]=1; break;
			case TS:		typeN[l]=2; break;
			case TY:		typeN[l]=3; break;
			case EC:		typeN[l]=4; break;
			case OTHERS:	typeN[l]=5; break;
			default: throw new IllegalArgumentException("unsupported type: "+types[l]);
			}
		}
		
		return typeN;
	}
	
	
	static class Scatter{
		//
		private int vcount =24;
		private int tcount =0;
		private int pointer=0;
		
		private int[] counters=null;
		
		private float[] seps=null;
		
		private Variable[] vs=null;
		
		
		//
		public Scatter(int tlength,float[] seps){
			tcount=tlength;
			
			vs=new Variable[vcount];
			counters=new int[seps.length-1];
			this.seps=seps;
			
			for(int i=0;i<vcount;i++)
			vs[i]=new Variable("vs"+i,false,new Range(tcount,1,1,1));
			
			vs[ 0].setName("sst"  );	vs[ 0].setCommentAndUnit("sea surface temperature (K)");
			vs[ 1].setName("vws"  );	vs[ 1].setCommentAndUnit("vertical wind shear (m/s)");
			vs[ 2].setName("mpi"  );	vs[ 2].setCommentAndUnit("maximum potential intensity (m/s)");
			vs[ 3].setName("efcs" );	vs[ 3].setCommentAndUnit("300-600 km eddy angular momentum flux convergence (m/s/day)");
			vs[ 4].setName("efcl" );	vs[ 4].setCommentAndUnit("500-800 km eddy angular momentum flux convergence (m/s/day)");
			vs[ 5].setName("wind" );	vs[ 5].setCommentAndUnit("10-min maximum sustained wind speed");
			vs[ 6].setName("rwnd" );	vs[ 6].setCommentAndUnit("storm-relative maximum sustained wind speed");
			vs[ 7].setName("pres" );	vs[ 7].setCommentAndUnit("minimum sea level pressure");
			vs[ 8].setName("dsst" );	vs[ 8].setCommentAndUnit("delta SST");
			vs[ 9].setName("dvws" );	vs[ 9].setCommentAndUnit("delta VWS");
			vs[10].setName("dmpi" );	vs[10].setCommentAndUnit("delta MPI");
			vs[11].setName("defcs");	vs[11].setCommentAndUnit("delta EFCs");
			vs[12].setName("defcl");	vs[12].setCommentAndUnit("delta EFCl");
			vs[13].setName("dwind");	vs[13].setCommentAndUnit("delta wind");
			vs[14].setName("drwnd");	vs[14].setCommentAndUnit("delta relative wind");
			vs[15].setName("dpres");	vs[15].setCommentAndUnit("delta pressure");
			vs[16].setName("uvel" );	vs[16].setCommentAndUnit("x-component of moving velocity");
			vs[17].setName("vvel" );	vs[17].setCommentAndUnit("y-component of moving velocity");
			vs[18].setName("mspd" );	vs[18].setCommentAndUnit("magnitude of moving velocity");
			vs[19].setName("lons" );	vs[19].setCommentAndUnit("longitudes");
			vs[20].setName("lats" );	vs[20].setCommentAndUnit("latitudes");
			vs[21].setName("magm" );	vs[21].setCommentAndUnit("magnitude of upper-level wind speed within 500 km");
			vs[22].setName("type" );	vs[22].setCommentAndUnit("type. 1-5 for TD, TS, TY, EC and others");
			vs[23].setName("isbm" );	vs[23].setCommentAndUnit("mean inertial stability");
		}
		
		
		// getor and setor
		public int getCount(){ return vcount;}
		
		public int getPointer(){ return pointer;}
		
		public Variable[] getVariables(){ return vs;}
		
		public void setPointer(int p){ pointer=p;}
		
		
		//
		public void compile(boolean[] valid,float[]... data){
			if(data.length!=vcount) throw new IllegalArgumentException("13 sets of data required");
			
			int N=data[0].length;
			
			if(N!=valid.length) throw new IllegalArgumentException("invalid data length");
			
			for(int i=1,I=data.length;i<I;i++)
			if(N!=data[i].length)
				throw new IllegalArgumentException("invalid data length "+N+" v.s. "+data[i].length+"\t"+i);
			
			for(int l=0;l<N;l++){
				if(valid[l]) for(int i=0;i<vcount;i++) vs[i].getData()[0][0][0][pointer]=data[i][l];
				else         for(int i=0;i<vcount;i++) vs[i].getData()[0][0][0][pointer]=undef;
				
				pointer++;
			}
		}
		
		public void count(){
			/*** count delta-pressure ***/
			if(!vs[15].getName().equalsIgnoreCase("dpres")) throw new IllegalArgumentException("incorrect variable");
			
			float[] dprs=vs[15].getData()[0][0][0];
			
			int validC=0;
			
			for(int l=0,L=dprs.length;l<L;l++)
			if(dprs[l]!=undef){
				validC++; boolean hasCount=false;
				
				for(int m=1,M=seps.length;m<M;m++) if(dprs[l]<seps[m]&&dprs[l]>=seps[m-1]){ counters[m-1]++; hasCount=true; break;}
				
				if(!hasCount) throw new IllegalArgumentException(dprs[l]+" is not counted "+(dprs[l]<seps[1])+"\t"+(dprs[l]>=seps[0])+"\t"+seps[0]);
			}
			
			System.out.println("delta pressure ("+deltaP_interval+" interval) stat: "+validC+" valid count");
			for(int m=1,M=seps.length;m<M;m++)
			System.out.println("  within ["+seps[m-1]+" - "+seps[m]+"): "+counters[m-1]);
		}
		
		//
		public void writeToFile(String path){
			CtlDataWriteStream cdws=new CtlDataWriteStream(path);
			cdws.writeData(vs);	cdws.closeFile();
			
			StringBuilder sb=new StringBuilder();
			
			sb.append("dset ^Scatter.dat\n");
			sb.append("undef "+undef+"\n");
			sb.append("title Best Track Scatter\n");
			sb.append("xdef     1 linear 0 0.02\n");
			sb.append("ydef     1 linear 0 0.02\n");
			sb.append("zdef     1 levels 1000\n");
			sb.append("tdef "+tcount+" linear 00:00Z01JAN1989 6hr\n");
			sb.append("vars "+vcount+"\n");
			for(int i=0;i<vcount;i++)
			sb.append(String.format("%-6s",vs[i].getName())+" 0 99 "+vs[i].getCommentAndUnit()+"\n");
			sb.append("endvars\n");
			
			try{
				FileWriter fw=new FileWriter(path.replace(".dat",".ctl"));
				fw.write(sb.toString());	fw.close();
				
			}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		}
	}
}
