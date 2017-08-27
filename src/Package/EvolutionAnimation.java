package Package;
//
import java.io.FileWriter;
import java.io.IOException;
import miniufo.basic.InterpolationModel.Type;
import miniufo.database.AccessBestTrack;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.MDate;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;

//
public class EvolutionAnimation{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Elena/Data.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		Range r=new Range("lon(240,300);lat(10,55);t(5,29)",dd);
		Variable[] wind=df.getVariables(r,"u","v");
		
		int t=wind[0].getTCount();
		
		Variable u=wind[0].interpolateT((t-1)*3+1,Type.LINEAR);	u.setUndef(-9999);
		Variable v=wind[1].interpolateT((t-1)*3+1,Type.LINEAR);	v.setUndef(-9999);
		
		Typhoon tr=AccessBestTrack.getTyphoonsFromNHC(
			"D:/Data/Typhoons/NHC/original/tracks1851to2008_atl_reanal.txt","name=Elena;time=1Aug1985-31Sep1985"
		).get(0);
		
		tr=tr.interpolateAlongT(2);System.out.println(tr);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/DiagnosisVortex/Elena/animation.dat");
		dw.writeData(dd,u,v);	dw.closeFile();
		
		//writeGS(dd,tr,"d:/Data/DiagnosisVortex/Elena/animation.gs");
		
		Variable[] inten=AccessBestTrack.toIntensityVariables(tr);
		dw=DataIOFactory.getDataWrite(dd,"d:/Data/DiagnosisVortex/Elena/intensity.dat");
		dw.writeData(inten);	dw.closeFile();
		
		writeCtl(tr,"d:/Data/DiagnosisVortex/Elena/intensity.ctl");
	}
	
	public static void writeGS(DataDescriptor dd,Typhoon tr,String path){
		StringBuilder sb=new StringBuilder();
		
		sb.append("'open "+path.replace("gs","ctl")+"'\n");
		sb.append("'enable print "+path.replace("gs","gmf")+"'\n\n");
		
		sb.append("lons=\"");
		for(int i=0;i<tr.getTCount();i++) sb.append(tr.getXPositions()[i]+" ");
		sb.append("\"\n");
		
		sb.append("lats=\"");
		for(int i=0;i<tr.getTCount();i++) sb.append(tr.getYPositions()[i]+" ");
		sb.append("\"\n\n");
		
		sb.append("'set grads off'\n");
		sb.append("'set grid off'\n");
		sb.append("'set lon 90 135'\n");
		sb.append("'set lat 15 60'\n\n");
		
		sb.append("cc=1\n");
		sb.append("tt="+(dd.getTNum(tr.getTime(0))*3+1)+"\n");
		sb.append("while(tt<="+(dd.getTNum(tr.getTime(tr.getTCount()-20))*3+1)+")\n");
		sb.append("'set t 'tt\n");
		sb.append("'setlopts 9 0.23 20 10'\n");
		sb.append("'set rbrange 0 60'\n");
		sb.append("'set cthick 5'\n");
		sb.append("'set arrowhead 0.07'\n");
		sb.append("'d u;v;mag(u,v)'\n");
		sb.append("lon=subwrd(lons,cc)\n");
		sb.append("lat=subwrd(lats,cc)\n");
		sb.append("'pty ' lon' 'lat' 0.5 1 9'\n");
		sb.append("'drawtime'\n");
		sb.append("'print'\n");
		sb.append("'c'\n");
		sb.append("cc=cc+1\n");
		sb.append("tt=tt+1\n");
		sb.append("endwhile\n\n");
		
		sb.append("'disable print'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		try(FileWriter fw=new FileWriter(path)){
			fw.write(sb.toString());
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	public static void writeCtl(Typhoon tr,String path){
		StringBuilder sb=new StringBuilder();
		
		long[] times=tr.getTimes();
		
		MDate str=new MDate(times[0]);
		
		int dtmili=str.getDT(new MDate(times[1]));
		int dt=dtmili/3600;
		
		if(dtmili%3600!=0) throw new IllegalArgumentException("invalid time increment");
		
		sb.append("dset "+path.replace("ctl","dat")+"\n");
		sb.append("undef -9999\n");
		sb.append("title JMA intensity of "+tr.getName()+"\n");
		sb.append("xdef 1 linear 0 1\n");
		sb.append("ydef 1 linear 0 1\n");
		sb.append("zdef 1 linear 0 1\n");
		sb.append("tdef "+tr.getTCount()+" linear "+str.toGradsDate()+" "+dt+"hr\n");
		sb.append("vars 2\n");
		sb.append("pres 0 99 winds (hPa)\n");
		sb.append("wnds 0 99 winds (m s^-1)\n");
		sb.append("endvars\n");
		
		try{
			FileWriter fw=new FileWriter(path);
			fw.write(sb.toString());	fw.close();
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
}
