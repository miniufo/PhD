package Package;
//
import miniufo.application.statisticsModel.EOFApplication;
import miniufo.application.statisticsModel.EOFResult;
import miniufo.basic.ArrayUtil;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;

//
public class CompositeEOF{
	//
	public static void main(String[] args){
		String path="d:/Data/PhD/Climatology/JMA/Composite/";
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"ComEFCS_SW.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] wind=df.getVariables(new Range("lev(200,200)",dd),false,"ru","rv");
		
		EOFResult reU=EOFApplication.EOF(wind[0],10,10);
		EOFResult reV=EOFApplication.EOF(wind[1],10,10);
		
		DataWrite dw=null;
		dw=DataIOFactory.getDataWrite(dd,path+"modSW.dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,reU.getModes(),reV.getModes()));	dw.closeFile();
		dw=DataIOFactory.getDataWrite(dd,path+"timSW.dat");
		dw.writeData(dd,ArrayUtil.concatAll(Variable.class,reU.getTimes(),reV.getTimes()));	dw.closeFile();
	}
}
