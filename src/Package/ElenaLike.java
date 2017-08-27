package Package;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import miniufo.application.advanced.CoordinateTransformation;
import miniufo.application.basic.DynamicMethodsInCC;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.CsmDescriptor;
import miniufo.descriptor.CtlDescriptor;
import miniufo.diagnosis.CylindricalSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.diagnosis.Variable.Dimension;
import miniufo.io.CtlDataWriteStream;
import miniufo.lagrangian.Typhoon;

//
public class ElenaLike{
	//
	private static int deltaP_interval=1;		// 1 for 6hr, 2 for 12hr...
	
	private static final DataSets dsets=DataSets.JMA;
	private static final String respath="d:/Data/PhD/Climatology/"+dsets+"/ElenaLike/";
	
	static final TCRec YURI=new TCRec(
		120f,   4f, 180f,  44f,"9128","YURI","1991",
		new MDate[]{new MDate(19911129060000L),new MDate(19911129180000L)}
	);
	
	static final TCRec OMAR=new TCRec(
		100f,   4f, 160f,  44f,"9215","OMAR","1992",
		new MDate[]{new MDate(19920831060000L),new MDate(19920831180000L)}
	);
	
	static final TCRec WARD=new TCRec(
		140f,   12f, 200f,  52f,"9221","WARD","1992",
		new MDate[]{new MDate(19920929180000L),new MDate(19920930060000L)}
	);
	
	static final TCRec KEONI=new TCRec(
		140f,   12f, 200f,  52f,"9310","KEONI","1993",
		new MDate[]{new MDate(19930822120000L),new MDate(19930823000000L)}
	);
	
	static final TCRec VIOLET=new TCRec(
		115f,   12f, 175f,  52f,"9617","VIOLET","1996",
		new MDate[]{new MDate(19960920180000L),new MDate(19960921060000L)}
	);
	
	static final TCRec KIROGI=new TCRec(
		115f,   8f, 175f,  48f,"0003","KIROGI","2000",
		new MDate[]{new MDate(20000706060000L),new MDate(20000706180000L)}
	);
	
	static final TCRec KONGREY=new TCRec(
		130f,   18f, 190f,  58f,"0106","KONG-REY","2001",
		new MDate[]{new MDate(20010727060000L),new MDate(20010727180000L)}
	);
	
	static final TCRec NESAT=new TCRec(
		110f,   3f, 170f,  43f,"0504","NESAT","2005",
		new MDate[]{new MDate(20050605180000L),new MDate(20050606060000L)}
	);
	
	static final TCRec MAWAR=new TCRec(
		115f,  12f, 175f,  52f,"0511","MAWAR","2005",
		new MDate[]{new MDate(20050822060000L),new MDate(20050822180000L)}
	);
	
	static final TCRec RYAN=new TCRec(
		120f,  12f, 180f,  52f,"9217","RYAN","1992",
		new MDate[]{new MDate(19920901180000L),new MDate(19920902060000L)}
	);
	
	static final TCRec MARIA=new TCRec(
		115f,  13f, 175f,  53f,"0607","MARIA","2006",
		new MDate[]{new MDate(20060804000000L),new MDate(20060804120000L)}
	);
	
	static final TCRec SHANSHAN=new TCRec(
		100f,   4f, 160f,  44f,"0613","SHANSHAN","2006",
		new MDate[]{new MDate(20060909000000L),new MDate(20060919180000L)}
	);
	
	
	//
	public static void main(String[] args){
		//processOne(YURI);
		//processOne(OMAR);
		//processOne(WARD);
		//processOne(KEONI);
		//processOne(VIOLET);
		//processOne(KIROGI);
		//processOne(KONGREY);
		//processOne(NESAT);
		//processOne(MAWAR);
		//processOne(RYAN);
		//processOne(MARIA);
		processOne(SHANSHAN);
	}
	
	static void processOne(TCRec rec){
		List<Typhoon> ls=AccessBestTrack.getTyphoons(IntensityModel.getPath(dsets),"time=1Jan"+rec.year+"-31Dec"+rec.year+";name="+rec.name,dsets);
		
		if(ls.size()!=1){ System.out.println("there are more than 1 records found"); System.exit(0);}
		
		CtlDescriptor ctl=(CtlDescriptor)DiagnosisFactory.getDataDescriptor("D:/Data/ERAInterim/Data.ctl");
		
		Typhoon tr=ls.get(0);
		
		DiagnosisFactory df=DiagnosisFactory.parseContent(tr.toCSMString("d:/ctl.ctl",36,25,2,0.3f,-650,850));
		CsmDescriptor dd=(CsmDescriptor)df.getDataDescriptor();
		
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
		
		Variable utm=utvr[0].anomalizeX();	utvr[1].anomalizeX();
		Variable efcsm=dm.cREFC(utvr[0],utvr[1]).averageAlong(Dimension.Y, 9,18);	// 300-600 km
		Variable etam=dm.cMeanAbsoluteVorticity(utm).averageAlong(Dimension.Y,9,18);
		
		float[] EFCs=efcsm.getData()[1][0][0];
		float[] VWSm= vwsm.getData()[0][0][0];
		float[] SSTm= sstm.getData()[0][0][0];
		float[] ETAm= etam.getData()[1][0][0];
		float[] pres=tr.getPressures();
		float[] dprs=Typhoon.getChangesByForwardDiff(pres,deltaP_interval);
		float[] mpim=IntensityModel.cMPI(SSTm);
		
		Range r=new Range(tr.getTCount(),1,1,1);
		
		Variable sst=new Variable("sst",false,r);
		Variable vws=new Variable("vws",false,r);
		Variable efc=new Variable("efc",false,r);
		Variable eta=new Variable("eta",false,r);
		Variable prs=new Variable("prs",false,r);
		Variable del=new Variable("del",false,r);
		Variable mpi=new Variable("mpi",false,r);
		
		System.arraycopy(SSTm,0,sst.getData()[0][0][0],0,tr.getTCount());
		System.arraycopy(VWSm,0,vws.getData()[0][0][0],0,tr.getTCount());
		System.arraycopy(EFCs,0,efc.getData()[0][0][0],0,tr.getTCount());
		System.arraycopy(ETAm,0,eta.getData()[0][0][0],0,tr.getTCount());
		System.arraycopy(pres,0,prs.getData()[0][0][0],0,tr.getTCount());
		System.arraycopy(dprs,0,del.getData()[0][0][0],0,tr.getTCount());
		System.arraycopy(mpim,0,mpi.getData()[0][0][0],0,tr.getTCount());
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(respath+rec.ID+rec.name+".dat");
		cdws.writeData(sst,vws,efc,eta,prs,del,mpi); cdws.closeFile();
		
		
		/**** write ctl ****/
		StringBuilder sb=new StringBuilder();
		
		sb.append("dset ^"+rec.ID+rec.name+".dat\n");
		sb.append("title TC case\n");
		sb.append("undef "+ctl.getUndef(null)+"\n");
		sb.append("xdef 1 linear 0 1\n");
		sb.append("ydef 1 linear 0 1\n");
		sb.append("zdef 1 linear 0 1\n");
		sb.append("tdef "+tr.getTCount()+" linear "+new MDate(tr.getTime(0)).toGradsDate()+" 6hr\n");
		sb.append("vars 7\n");
		sb.append("sst 0 99 sea surface temperature within 500km\n");
		sb.append("vws 0 99 vertical wind shear within 500km\n");
		sb.append("efc 0 99 eddy momentum flux convergence within 300-600km\n");
		sb.append("eta 0 99 absolute vorticity within 300-600 km\n");
		sb.append("prs 0 99 minimum sea level pressure\n");
		sb.append("del 0 99 6-hourly delta-P\n");
		sb.append("mpi 0 99 maximum potential intensity determined by SST\n");
		sb.append("endvars\n");
		
		try(FileWriter fw=new FileWriter(respath+rec.ID+rec.name+".ctl")){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
		
		
		/**** write gs ****/
		StringBuilder gs=new StringBuilder();
		gs.append(rec.name+"lons=\""); for(int i=0;i<tr.getTCount();i++) gs.append(tr.getXPosition(i)+" "); gs.append("\"\n");
		gs.append(rec.name+"lats=\""); for(int i=0;i<tr.getTCount();i++) gs.append(tr.getYPosition(i) +" "); gs.append("\"\n");
		gs.append(rec.name+"type=\""); for(int i=0;i<tr.getTCount();i++) gs.append(tr.getTypes()[i]  +" "); gs.append("\"\n\n");
		
		gs.append("'open D:/Data/ERAInterim/Data.ctl'\n");
		gs.append("'enable print "+respath+rec.ID+rec.name+".gmf'\n\n");
		
		gs.append("'setvpage 4 3 4 1'\n");
		gs.append("'setlopts 5 0.22 15 10'\n");
		gs.append("'set grid off'\n");
		gs.append("'set time "+rec.tims[0].toGradsDate()+"'\n");
		gs.append("'set lon "+rec.lonmin+" "+rec.lonmax+"'\n");
		gs.append("'set lat "+rec.latmin+" "+rec.latmax+"'\n");
		gs.append("'set lev 200'\n");
		gs.append("'set arrowhead -0.35'\n");
		gs.append("'set ccolor 15'\n");
		gs.append("'set arrscl 0.4 40'\n");
		gs.append("'d maskout(u,40-mag(u,v));v'\n");
		gs.append("'set ccolor 1'\n");
		gs.append("'set arrscl 0.4 100'\n");
		gs.append("'d maskout(u,mag(u,v)-40);v'\n");
		gs.append("lon=subwrd("+rec.name+"lons,"+(tr.getTag(rec.tims[0].getLongTime())+1)+")\n");
		gs.append("lat=subwrd("+rec.name+"lats,"+(tr.getTag(rec.tims[0].getLongTime())+1)+")\n");
		gs.append("typ=subwrd("+rec.name+"type,"+(tr.getTag(rec.tims[0].getLongTime())+1)+")\n");
		gs.append("'set line 1 1 1'\n");
		gs.append("'drawmark 3 'lon' 'lat' 0.25'\n");
		gs.append("tt=1\n");
		gs.append("while(tt<"+tr.getTCount()+")\n");
		gs.append("lon1=subwrd("+rec.name+"lons,tt)\n");
		gs.append("lat1=subwrd("+rec.name+"lats,tt)\n");
		gs.append("typ1=subwrd("+rec.name+"type,tt)\n");
		gs.append("lon2=subwrd("+rec.name+"lons,tt+1)\n");
		gs.append("lat2=subwrd("+rec.name+"lats,tt+1)\n");
		gs.append("typ2=subwrd("+rec.name+"type,tt+1)\n");
		gs.append("'set line 1 1 2'\n");
		gs.append("'drawline 'lon1' 'lat1' 'lon2' 'lat2\n");
		gs.append("if(typ1='TD')\n"   ); gs.append("'set line 4 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TS')\n"   ); gs.append("'set line 3 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TY')\n"   ); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='EC')\n"   ); gs.append("'set line 7 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='OTHER')\n"); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("'drawmark 3 'lon1' 'lat1' 0.1'\n");
		gs.append("tt=tt+1\n");
		gs.append("endwhile\n");
		gs.append("if(typ1='TD')\n"   ); gs.append("'set line 4 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TS')\n"   ); gs.append("'set line 3 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TY')\n"   ); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='EC')\n"   ); gs.append("'set line 7 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='OTHER')\n"); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("'drawmark 3 'lon2' 'lat2' 0.1'\n");
		gs.append("'drawtime'\n\n");
		
		gs.append("'setvpage 4 3 4 2'\n");
		gs.append("'setlopts 5 0.22 15 10'\n");
		gs.append("'set time "+rec.tims[1].toGradsDate()+"'\n");
		gs.append("'set arrowhead -0.35'\n");
		gs.append("'set ccolor 15'\n");
		gs.append("'set arrscl 0.4 40'\n");
		gs.append("'d maskout(u,40-mag(u,v));v'\n");
		gs.append("'set ccolor 1'\n");
		gs.append("'set arrscl 0.4 100'\n");
		gs.append("'d maskout(u,mag(u,v)-40);v'\n");
		gs.append("lon=subwrd("+rec.name+"lons,"+(tr.getTag(rec.tims[1].getLongTime())+1)+")\n");
		gs.append("lat=subwrd("+rec.name+"lats,"+(tr.getTag(rec.tims[1].getLongTime())+1)+")\n");
		gs.append("typ=subwrd("+rec.name+"type,"+(tr.getTag(rec.tims[1].getLongTime())+1)+")\n");
		gs.append("'set line 1 1 1'\n");
		gs.append("'drawmark 3 'lon' 'lat' 0.25'\n");
		gs.append("tt=1\n");
		gs.append("while(tt<"+tr.getTCount()+")\n");
		gs.append("lon1=subwrd("+rec.name+"lons,tt)\n");
		gs.append("lat1=subwrd("+rec.name+"lats,tt)\n");
		gs.append("typ1=subwrd("+rec.name+"type,tt)\n");
		gs.append("lon2=subwrd("+rec.name+"lons,tt+1)\n");
		gs.append("lat2=subwrd("+rec.name+"lats,tt+1)\n");
		gs.append("typ2=subwrd("+rec.name+"type,tt+1)\n");
		gs.append("'set line 1 1 2'\n");
		gs.append("'drawline 'lon1' 'lat1' 'lon2' 'lat2\n");
		gs.append("if(typ1='TD')\n"   ); gs.append("'set line 4 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TS')\n"   ); gs.append("'set line 3 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TY')\n"   ); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='EC')\n"   ); gs.append("'set line 7 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='OTHER')\n"); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("'drawmark 3 'lon1' 'lat1' 0.1'\n");
		gs.append("tt=tt+1\n");
		gs.append("endwhile\n");
		gs.append("if(typ1='TD')\n"   ); gs.append("'set line 4 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TS')\n"   ); gs.append("'set line 3 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='TY')\n"   ); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='EC')\n"   ); gs.append("'set line 7 1 1'\n"); gs.append("endif\n");
		gs.append("if(typ1='OTHER')\n"); gs.append("'set line 2 1 1'\n"); gs.append("endif\n");
		gs.append("'drawmark 3 'lon2' 'lat2' 0.1'\n");
		gs.append("'drawtime'\n\n");
		gs.append("'close 1'\n\n");
		
		gs.append("'open "+respath+rec.ID+rec.name+".ctl'\n");
		gs.append("'setvpage 4 3.6 4 3.5'\n");
		gs.append("'setlopts 5 0.25 1 10'\n");
		gs.append("'set parea 1 8 0.87 6.73'\n");
		gs.append("'set grid on'\n");
		gs.append("'set t 1 last'\n");
		gs.append("'set vrange 880 980'\n");
		gs.append("'set ccolor 1 1 4'\n");
		gs.append("'set cmark 5'\n");
		gs.append("'set digsize 0.14'\n");
		gs.append("'d prs'\n");
		gs.append("'set vrange 25 30'\n");
		gs.append("'set ylpos 0 r'\n");
		gs.append("'set ylint 1'\n");
		gs.append("'set ccolor 2 1 4'\n");
		gs.append("'set cmark 0'\n");
		gs.append("'d sst'\n");
		gs.append("'set vrange 0 15'\n");
		gs.append("'set ylpos 0.7 r'\n");
		gs.append("'set ylint 3'\n");
		gs.append("'set ccolor 4 1 4'\n");
		gs.append("'set cmark 0'\n");
		gs.append("'d vws'\n");
		gs.append("'set vrange -10 40'\n");
		gs.append("'set ylpos 1.4 r'\n");
		gs.append("'set ylint 10'\n");
		gs.append("'set ccolor 3 1 4'\n");
		gs.append("'set cmark 0'\n");
		gs.append("'d efc'\n");
		gs.append("'draw title "+rec.name+" ("+rec.year+") SST VWS EFC'\n\n");
		
		gs.append("'print'\n");
		gs.append("'c'\n\n");
		
		gs.append("'close 1'\n");
		gs.append("'disable print'\n");
		gs.append("'reinit'\n");
		
		try(FileWriter fw=new FileWriter(respath+rec.ID+rec.name+".gs")){
			fw.write(gs.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	
	static final class TCRec{
		//
		float lonmin=0;
		float latmin=0;
		float lonmax=0;
		float latmax=0;
		
		String ID=null;
		String name=null;
		String year=null;
		
		MDate[] tims=null;
		
		public TCRec
		(float lonmin,float latmin,float lonmax,float latmax,String ID,String name,String year,MDate[] tims){
			this.lonmin=lonmin;
			this.latmin=latmin;
			this.lonmax=lonmax;
			this.latmax=latmax;
			this.ID=ID;
			this.name=name;
			this.year=year;
			this.tims=tims;
		}
	}
}
