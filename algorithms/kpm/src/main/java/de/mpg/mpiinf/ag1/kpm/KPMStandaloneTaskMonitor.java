package de.mpg.mpiinf.ag1.kpm;

import dk.sdu.kpm.taskmonitors.IKPMTaskMonitor;

public class KPMStandaloneTaskMonitor implements IKPMTaskMonitor {
	private int nextToPrint = 0;

	@Override
	public void setTitle(String title) {
        System.out.println(title);}

	@Override
	public void setProgress(double progress) {
		progress *= 100;

		//if(progress > nextToPrint){
			StringBuilder bar = new StringBuilder("[");

			for(int i = 0; i < 50; i++){
				if( i < (progress/2)){
					bar.append("=");
				}else if( i == (progress/2)){
					bar.append(">");
				}else{
					bar.append(" ");
				}
			}

			bar.append("]   " + progress + "%\r");
			System.out.print(bar.toString());
			//nextToPrint += 10;
		//}
	}

	@Override
	public void setStatusMessage(String statusMessage) {
		System.out.println(statusMessage);
	}

}
