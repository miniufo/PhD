package Package;
//
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import miniufo.descriptor.CsmDescriptor;
import miniufo.diagnosis.DiagnosisFactory;

//
public class CylindricalGrids{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/Haima/Haima.csm");
		CsmDescriptor csd=(CsmDescriptor)df.getDataDescriptor();
		
		float[][][] lons=csd.getLon();
		float[][][] lats=csd.getLat();
		
		StringBuilder sb=new StringBuilder();
		
		sb.append("'open d:/Data/DiagnosisVortex/Haima/Haima.ctl'\n");
		sb.append("'enable print d:/Paper/PhD/Chapter2/Data/CylindricalGrids.gmf'\n\n");
		
		int tt=19;
		sb.append("'setvpage 2 2 2 1'\n");
		sb.append("'setlopts 8 0.26 10 10'\n");
		sb.append("'set mpdset mres'\n");
		sb.append("'set map 1 1 1'\n");
		sb.append("'set lon 93 140'\n");
		sb.append("'set lat 28 63'\n");
		sb.append("'set t "+tt+"'\n");
		sb.append("'set cmin 10000'\n");
		sb.append("'d t'\n");
		sb.append("'drawtime'\n");
		for(int j=1;j<csd.getYCount();j+=2)
		for(int i=0;i<csd.getXCount();i+=2)
		sb.append("'drawmark 3 "+lons[tt-1][j][i]+" "+lats[tt-1][j][i]+" 0.08'\n");
		sb.append("'drawmark 3 "+lons[tt-1][csd.getYCount()-1][0]+" "+lats[tt-1][csd.getYCount()-1][0]+" 0.08'\n\n");
		
		sb.append("'setvpage 2 2 2 2'\n");
		sb.append("'setlopts 8 0.26 10 10'\n");
		sb.append("'set map 1 1 1'\n");
		sb.append("'set mpdset mres'\n");
		sb.append("'set mproj lambert'\n");
		sb.append("'set lon 90 143'\n");
		sb.append("'set lat 28 63'\n");
		sb.append("'set t "+tt+"'\n");
		sb.append("'set cmin 10000'\n");
		sb.append("'d t'\n");
		sb.append("'drawtime'\n");
		for(int j=1;j<csd.getYCount();j+=2)
		for(int i=0;i<csd.getXCount();i+=2)
		sb.append("'drawmark 3 "+lons[tt-1][j][i]+" "+lats[tt-1][j][i]+" 0.08'\n");
		sb.append("'drawmark 3 "+lons[tt-1][csd.getYCount()-1][0]+" "+lats[tt-1][csd.getYCount()-1][0]+" 0.08'\n\n");
		
		sb.append("'print'\n");
		sb.append("'c'\n\n");
		sb.append("'disable print'\n");
		sb.append("'close 1'\n");
		sb.append("'reinit'\n");
		
		toFile(sb,"d:/Paper/PhD/Chapter2/Data/CylindricalGrids.gs");
	}
	
	static void toFile(StringBuilder sb,String path){
		try{
			FileWriter fw=new FileWriter(new File(path));
			fw.write(sb.toString());	fw.close();
			
		}catch(IOException e){e.printStackTrace(); System.exit(0);}
	}
}
