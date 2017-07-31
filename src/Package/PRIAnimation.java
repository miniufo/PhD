//
package Package;

import java.io.File;
import java.io.FileWriter;
import miniufo.application.basic.IndexInSC;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel.Type;
import miniufo.database.AccessBestTrack;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;

//
public class PRIAnimation{
	//
	private static final String name="Linfa";
	private static final String year="2003";
	private static final String rngs="lon(105,165);lat(9,61)";
	private static final String path="d:/Data/DiagnosisVortex/";
	
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+name+"/"+name+".ctl");
		DataDescriptor dd=df.getDataDescriptor();
		Range r=new Range("lev(200,200)",dd);
		Variable[] wind=df.getVariables(r,"u","v");
		
		int t=wind[0].getTCount();
		
		Variable u=wind[0].interpolateT((t-1)*3+1,Type.LINEAR);	u.setUndef(-9999);
		Variable v=wind[1].interpolateT((t-1)*3+1,Type.LINEAR);	v.setUndef(-9999);
		
		DataWrite dw1=DataIOFactory.getDataWrite(dd,path+name+"/PRI/Data.dat");
		dw1.writeData(dd,u,v);	dw1.closeFile();
		
		
		DiagnosisFactory df2=DiagnosisFactory.parseFile(path+name+"/PRI/Data.ctl");
		DataDescriptor dd2=df2.getDataDescriptor();
		Range r2=new Range(rngs,dd2);
		Variable[] wind2=df2.getVariables(r2,"u","v");
		
		Typhoon tr=AccessBestTrack.getTyphoonsFromJMA(
			"D:/Data/Typhoons/JMA/JMA.txt","name="+name+";time=1Jan"+year+"-31Dec"+year
		).get(0);
		
		System.out.println(tr);
		
		tr=tr.interpolateAlongT(2);
		
		Range rBST=new Range(tr.getTCount(),1,1,1);
		
		Variable wnd=new Variable("wnd",false,rBST);	wnd.setUndef(-9999);
		Variable prs=new Variable("prs",false,rBST);	prs.setUndef(-9999);
		
		float[][][][] wdata=wnd.getData();
		float[][][][] pdata=prs.getData();
		
		System.arraycopy(tr.getWinds(),0,wdata[0][0][0],0,tr.getTCount());
		System.arraycopy(tr.getPressures(),0,pdata[0][0][0],0,tr.getTCount());
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+name+"/PRI/intensity.dat");
		dw.writeData(dd,wnd,prs);	dw.closeFile();
		
		Variable[] idx1=IndexInSC.c2DHorizontalIndex(
			dd2,rngs,tr,0.3f,19,36,"REFC"//,"PEFC","AEFC","EAMA","FFCT","FFBS"
		);
		//Variable[] idx2=IndexInSC.c2DHorizontalIndex(
		//	dd2,"lon(70,180);lat(12,58)",0.3f,19,72,"REFC","PEFC","AEFC","EAMA","FFCT","FFBS"
		//);
		
		DataWrite dw2=DataIOFactory.getDataWrite(dd2,path+name+"/PRI/PRIanimation.dat");
		dw2.writeData(dd2,ArrayUtil.concatAll(Variable.class,wind2,idx1));	dw2.closeFile();
		writeGS(path,tr,dd2);
	}
	
	public static void writeGS(String path,Typhoon tr,DataDescriptor dd){
		String[] lonlat=rngs.split(";");
		String lons=lonlat[0].replace("lon(","").replace(")","").replace(","," ");
		String lats=lonlat[1].replace("lat(","").replace(")","").replace(","," ");
		
		StringBuffer sb=new StringBuffer();
		
		sb.append("'open "+path+name+"/PRI/PRIanimation.ctl'\n");
		sb.append("'open "+path+name+"/PRI/intensity.ctl'\n\n");
		sb.append("lons=\""); for(int l=0;l<tr.getTCount();l++) sb.append(tr.getLongitudes()[l]+" "); sb.append("\"\n");
		sb.append("lats=\""); for(int l=0;l<tr.getTCount();l++) sb.append(tr.getLongitudes()[l]+" "); sb.append("\"\n\n");
		sb.append("'set rgb 16   0   0 255'\n");
		sb.append("'set rgb 17  55  55 255'\n");
		sb.append("'set rgb 18 110 110 255'\n");
		sb.append("'set rgb 19 165 165 255'\n");
		sb.append("'set rgb 20 220 220 255'\n\n");
		sb.append("'set rgb 21 255 220 220'\n");
		sb.append("'set rgb 22 255 165 165'\n");
		sb.append("'set rgb 23 255 110 110'\n");
		sb.append("'set rgb 24 255  55  55'\n");
		sb.append("'set rgb 25 255   0   0'\n\n");
		sb.append("'set grads off'\n");
		sb.append("'set mpdset mres'\n");
		sb.append("'set map 1 1 1'\n\n");
		sb.append("cc=1\n");
		sb.append("tt=1\n");
		sb.append("while(tt<="+dd.getTCount()+")\n");
		sb.append("'set t 'tt\n");
		sb.append("'q time'\n");
		sb.append("tim=subwrd(result,3)\n");
		sb.append("'set lon "+lons+"'\n");
		sb.append("'set lat "+lats+"'\n");
		sb.append("'set lev 200'\n");
		sb.append("'setvpage 1.4 2 1.4 1.04'\n");
		sb.append("'setlopts 9 0.17 10 10'\n");
		sb.append("'set grid off'\n");
		sb.append("'set gxout shaded'\n");
		sb.append("'set clevs -80 -60 -40 -20 20 40 60 80'\n");
		sb.append("'set ccols 16 17 18 19 0 22 23 24 25'\n");
		sb.append("'d AEFC'\n");
		sb.append("'set gxout contour'\n");
		sb.append("'set rbrange 5 55'\n");
		sb.append("'set cthick 6'\n");
		sb.append("'set arrowhead -0.3'\n");
		sb.append("'set arrscl 0.6 70'\n");
		sb.append("'set arrlab off'\n");
		sb.append("'d skip(u,2);v'\n");
		sb.append("'cbarn 1 0 3.9 0.18'\n");
		sb.append("lon=subwrd(lons,cc)\n");
		sb.append("lat=subwrd(lats,cc)\n");
		sb.append("'pty ' lon' 'lat' 0.4 1 9'\n");
		sb.append("'drawtime'\n");
		sb.append("'set lon 0'\n");
		sb.append("'set lat -90'\n");
		sb.append("'set lev 1000'\n");
		sb.append("'d prs.2'\n");
		sb.append("prs=subwrd(result,4)\n\n");
		sb.append("'setvpage 1.9 2 1.578 1.87'\n");
		sb.append("'setlopts 9 0.22 1 4'\n");
		sb.append("'set grid on'\n");
		sb.append("'set lon 0'\n");
		sb.append("'set lat -90'\n");
		sb.append("'set lev 1000'\n");
		sb.append("'set cthick 9'\n");
		sb.append("'set t 1 "+dd.getTCount()+"'\n");
		sb.append("'set cmark 0'\n");
		sb.append("'d prs.2'\n");
		sb.append("'q w2xy 'tim' 'prs\n");
		sb.append("xx=subwrd(result,3)\n");
		sb.append("yy=subwrd(result,6)\n");
		sb.append("'draw mark 3 'xx' 'yy' 0.3'\n");
		sb.append("'draw title minimum central pressure'\n");
		sb.append("ttt=tt+100\n");
		sb.append("'printim "+path+name+"/PRI/PRIanimation'ttt'.png white x1024 y768'\n");
		sb.append("'c'\n");
		sb.append("cc=cc+1\n");
		sb.append("tt=tt+1\n");
		sb.append("endwhile\n\n");
		sb.append("'disable print'\n");
		sb.append("'close 2'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n\n");
		sb.append("'!convert -delay 10 -crop 930x510+35+55 +repage "+path+name+"/PRI/PRIanimation*.png "+path+name+"/PRI/PRIanimation.gif'\n");
		sb.append("'!rm "+path+name+"/PRI/PRIanimation*.png'\n");
		
		try{
			FileWriter fw=new FileWriter(new File(path+name+"/PRI/PRIanimation.gs"));
			fw.write(sb.toString());	fw.close();
			
		}catch(Exception e){ e.printStackTrace(); System.exit(0);}
	}
}
