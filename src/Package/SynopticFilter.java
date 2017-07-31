package Package;
//
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataRead;
import miniufo.io.DataWrite;
import miniufo.statistics.FilterModel;

//
public class SynopticFilter{
	//
	public static void main(String[] args){
		DataDescriptor ctl=DiagnosisFactory.getDataDescriptor("d:/Filter.ctl");
		
		Range r=new Range("y(1,1)",ctl);
		
		Variable v=new Variable("u",r);
		
		float[] buffer=new float[ctl.getTCount()];
		float[][] vdata=v.getData()[0][0];
		
		DataRead  dr=DataIOFactory.getDataRead(ctl);
		DataWrite dw=DataIOFactory.getDataWrite(ctl,"d:/SFilterd.dat");
		
		for(int j=0;j<ctl.getYCount();j++){
			r.setYRange(j+1);
			System.out.println("reading "+j);
			dr.readData(v);
			System.out.println("finish");
			for(int i=0;i<ctl.getXCount();i++){
				FilterModel.ButterworthFilter(vdata[i],buffer,8,32);
				System.arraycopy(buffer,0,vdata[i],0,ctl.getTCount());
			}
			
			dw.writeData(v);
		}
		
		dw.closeFile();	dr.closeFile();
	}
}
