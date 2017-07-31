//
package Package;

import miniufo.application.basic.IndexInSC;
import miniufo.basic.ArrayUtil;
import miniufo.database.AccessBestTrack;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.lagrangian.Typhoon;


//
public class HorizontalIndex{
	//
	private static final String dset="D:/Data/Typhoons/JMA/JMA.txt";
	private static final String name="Haima";
	private static final String cond="name="+name+";time=1Sep2004-30Sep2004";
	
	//
	public static void main(String[] args){
		//
		DiagnosisFactory df=DiagnosisFactory.parseFile("d:/Data/DiagnosisVortex/"+name+"/PRI/output.nc");
		DataDescriptor dd=df.getDataDescriptor();
		
		Typhoon tr=AccessBestTrack.getTyphoonsFromJMA(dset,cond).get(0);
		
		Range r=new Range(tr.getTRange(),dd);
		Variable[] wind=df.getVariables(r,"u","v");
		
		Variable[] idx1=IndexInSC.c2DHorizontalIndex(dd,"lon(90,150);lat(15,55);"+tr.getTRange(),tr,
		0.3f,19,72,"REFC","PEFC","AEFC","EAMA","FFCT","FFBS");
		Variable[] idx2=IndexInSC.c2DHorizontalIndex(dd,"lon(90,150);lat(15,55);"+tr.getTRange(),
		0.3f,19,72,"REFC","PEFC","AEFC","EAMA","FFCT","FFBS");
		
		DataWrite dw2=DataIOFactory.getDataWrite(dd,"d:/Data/DiagnosisVortex/"+name+"/PRI/PRIanimation.dat");
		dw2.writeData(dd,ArrayUtil.concatAll(Variable.class,wind,idx1,idx2));	dw2.closeFile();
	}
}
