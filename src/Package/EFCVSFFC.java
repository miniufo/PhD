package Package;
//
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
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;

//
public class EFCVSFFC{
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
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Climatology/"+dsets+"/Scatter/";
	private static final String tranges="time=1Jan1989-31Dec2009";
	
	//
	public static void main(String[] args){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),tranges,dsets);
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		
		Scatter sctr=new Scatter(18739);	// JMA total:18739, filtered:8642
		
		for(Typhoon tr:ls){
			DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString(36,28,2,0.3f,-650));
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
			
			Variable[] vars=df.getVariables(new Range("",dd),false,"u","v");
			
			Variable[] utvr=ct.reprojectToCylindrical(vars[0],vars[1]);
			dm.cStormRelativeAziRadVelocity(tr.getUVel(),tr.getVVel(),utvr[0],utvr[1]);
			
			utvr[0].anomalizeX();	utvr[1].anomalizeX();
			Variable EFCS=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),9 ,18);	// 300-600 km
			Variable EFCL=dm.cRadialAverage(dm.cREFC(utvr[0],utvr[1]),18,27);	// 600-900 km
			
			effectiveCount+=EFCS.getTCount();
			effectiveTC++;
			if(effectiveTC%20==0) System.out.print(".");
			
			sctr.compile(
				EFCS.getData()[1][0][0],	// EFCS
				EFCL.getData()[1][0][0]		// EFCL
			);
		}
		
		System.out.println("\n\npointer: "+sctr.getPointer());
		
		System.out.println();
		System.out.println("excluding depression : "+noDepress);
		System.out.println("excluding landing    : "+noLanding);
		System.out.println("excluding outside WNP: "+withinWNP);
		System.out.println("effective TCs  :"+effectiveTC);
		System.out.println("effective count:"+effectiveCount);
		
		sctr.writeToFile(respath+"ScatterEFCVSFFC.dat");
	}
	
	
	static class Scatter{
		//
		private int vcount =4;
		private int tcount =0;
		private int pointer=0;
		
		private Variable[] vs=null;
		
		
		//
		public Scatter(int tlength){
			tcount=tlength;
			
			vs=new Variable[vcount];
			
			for(int i=0;i<vcount;i++)
			vs[i]=new Variable("vs"+i,false,new Range(tcount,1,1,1));
			
			vs[0].setName("EFCS");	vs[0].setCommentAndUnit("EFCS");
			vs[1].setName("FFCS");	vs[1].setCommentAndUnit("FFCS");
			vs[2].setName("EFCL");	vs[2].setCommentAndUnit("EFCL");
			vs[3].setName("FFCL");	vs[3].setCommentAndUnit("FFCL");
		}
		
		
		// getor and setor
		public int getCount(){ return vcount;}
		
		public int getPointer(){ return pointer;}
		
		public Variable[] getVariables(){ return vs;}
		
		public void setPointer(int p){ pointer=p;}
		
		
		//
		public void compile(float[]... data){
			if(data.length!=vcount) throw new IllegalArgumentException("12 sets of data required");
			
			int N=data[0].length;
			
			for(int i=1,I=data.length;i<I;i++)
			if(N!=data[i].length)
			throw new IllegalArgumentException("invalid data length "+N+" v.s. "+data[i].length+"\t"+i);
			
			for(int l=0;l<N;l++){
				for(int i=0;i<vcount;i++)
				vs[i].getData()[0][0][0][pointer]=data[i][l];
				
				pointer++;
			}
		}
		
		//
		public void writeToFile(String path){
			CtlDataWriteStream cdws=new CtlDataWriteStream(path);
			cdws.writeData(vs);	cdws.closeFile();
			
			StringBuilder sb=new StringBuilder();
			
			sb.append("dset ^ScatterEFCVSFFC.dat\n");
			sb.append("undef -32767\n");
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
