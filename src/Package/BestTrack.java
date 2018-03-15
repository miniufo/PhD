//
package Package;

import java.util.List;
import miniufo.application.statisticsModel.BinningStatistics;
import miniufo.database.AccessBestTrack;
import miniufo.database.AccessBestTrack.DataSets;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;
import miniufo.lagrangian.Typhoon.TYPE;

//
public class BestTrack{
	//
	private static final int str=1989;
	private static final int end=2009;
	
	private static final DataSets[] ds={DataSets.CMA,DataSets.JMA,DataSets.JTWC};
	
	//
	public static void main(String[] args){
		DataDescriptor dd=DiagnosisFactory.getDataDescriptor("D:/Data/PhD/Climatology/Template.ctl");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"D:/Data/PhD/Climatology/BestTrack.dat");
		
		BinningStatistics bs=new BinningStatistics(dd);
		
		for(DataSets d:ds){
			List<Typhoon> ls=
			AccessBestTrack.getTyphoons("D:/Data/Typhoons/"+d+"/"+d+".txt","time=01Jan"+str+"-31Dec"+end,d);
			
			AccessBestTrack.recordsToFile(ls,"d:/Data/PhD/Climatology/"+d+".txt");
			
			Variable[] re1=bs.binningData(ls,Typhoon.UVEL,Typhoon.VVEL);
			Variable ace=bs.binningTCACE(ls,r->true);
			Variable foc=bs.binningCount(ls);
			Variable gfr=bs.binningTCGenesisFrequency(ls,17.2f);
			
			re1[0].setName("uts"+d);	dw.writeData(re1[0]);
			re1[1].setName("vts"+d);	dw.writeData(re1[1]);
			ace.setName("ace"+d);		dw.writeData(ace);
			foc.setName("foc"+d);		dw.writeData(foc);
			gfr.setName("gpt"+d);		dw.writeData(gfr);
			
			if(d!=DataSets.JTWC){
				Variable[] re2=bs.binningTCTypeCount(ls,r->true,TYPE.TD,TYPE.TS,TYPE.TY,TYPE.EC,TYPE.OTHERS);
				
				re2[0].setName("td"+d);	dw.writeData(re2[0]);
				re2[1].setName("ts"+d);	dw.writeData(re2[1]);
				re2[2].setName("ty"+d);	dw.writeData(re2[2]);
				re2[3].setName("ec"+d);	dw.writeData(re2[3]);
				re2[4].setName("ot"+d);	dw.writeData(re2[4]);
			}
		}
		
		dw.writeCtl(dd);	dw.closeFile();
	}
}
