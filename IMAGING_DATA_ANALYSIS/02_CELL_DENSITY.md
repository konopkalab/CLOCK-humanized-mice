# CELL DENSITY
- [2_CellDensity_FIJI.ijm](2_CellDensity_FIJI.ijm)

```
//open image and adjust its size and location
open()
name=getTitle;
run("Duplicate...", "title=dup.czi duplicate");
run("Set... ", "zoom=50");
setLocation(13, 13, 999, 999);
run("Set... ", "zoom=50");
setTool("hand");
waitForUser("Please go to ROI. Click OK when done");

//split and rename channels
x=1;
while (x==1) {
	setTool("rectangle");
	waitForUser("Please draw a rectangle to cover the window. Click OK when done");
	if (selectionType() == 0) {
		selectWindow("dup.czi");
		run("Duplicate...", "title=NeuN duplicate channels=1");
		selectWindow("dup.czi");
		run("Duplicate...", "title=Olig2 duplicate channels=2");
		selectWindow("dup.czi");
		run("Duplicate...", "title=CLOCK duplicate channels=3");
		selectWindow("dup.czi");
		run("Duplicate...", "title=DAPI duplicate channels=4");
		selectWindow("dup.czi");
		Overlay.addSelection("yellow", 19); 
		run("Make Composite");
		run("Stack to RGB");
		dir = getDirectory("choose a directory");
		saveAs("Tiff",dir+name+"_oPFC");
		close();
		x=0;
	}
}

//duplicate CLOCK channel for future usage
selectWindow("CLOCK");
run("Duplicate...", "title=CLOCK+Oligo");
run("RGB Color");
setForegroundColor(255, 0, 0);
run("Duplicate...", "title=CLOCK+Neuron");
run("RGB Color");
setForegroundColor(255, 0, 0);
run("Duplicate...", "title=CLOCK+Cell");
run("RGB Color");
setForegroundColor(255, 0, 0);

//generate blank images for data analysis
selectWindow("DAPI");
run("Duplicate...", "title=blackboard");
run("Make Binary");
selectWindow("blackboard");
run("Multiply...", "value=0.000");
run("Duplicate...", "title=DAPI_CLOCK");
run("Duplicate...", "title=Olig2_DAPI");
run("Duplicate...", "title=Olig2_DAPI_CLOCK");
run("Duplicate...", "title=NeuN_DAPI");
run("Duplicate...", "title=NueN_DAPI_CLOCK");

//DAPI image proccessing and mask generation
selectWindow("DAPI");
setLocation(1444, 13, 999, 999);
run("Set... ", "zoom=50");
run("Brightness/Contrast...");
waitForUser("Ajust brightest. Click OK when done");
run("Threshold...");
waitForUser("Ajust threshold. Click OK when done");
run("Subtract Background...");
run("Remove Outliers...");
run("Make Binary");
run("Options...");
run("Options...");
run("Duplicate...", "title=DAPI_mask");
selectWindow("DAPI_mask");
setLocation(1444, 13, 999, 999);
run("Set... ", "zoom=50");
run("Invert");

//CLOCK image proccessing and mask generation
selectWindow("CLOCK");
setLocation(1444, 13, 999, 999);
run("Set... ", "zoom=50");
run("Brightness/Contrast...");
waitForUser("Ajust brightest. Click OK when done");
run("Threshold...");
waitForUser("Ajust threshold. Click OK when done");
run("Subtract Background...");
run("Remove Outliers...");
run("Make Binary");
run("Options...");
run("Options...");
run("Duplicate...", "title=CLOCK_mask");
selectWindow("CLOCK_mask");
setLocation(1444, 13, 999, 999);
run("Set... ", "zoom=50");
run("Invert");

//Olig2 image proccessing
selectWindow("Olig2");
run("Set... ", "zoom=50");
setLocation(1444, 13, 999, 999);
run("Set... ", "zoom=50");
run("Brightness/Contrast...");
waitForUser("Ajust brightest. Click OK when done");
run("Threshold...");
waitForUser("Ajust threshold. Click OK when done");
run("Subtract Background...");
run("Remove Outliers...");
run("Make Binary");
run("Options...");
run("Options...");

//NeuN image proccessing
selectWindow("NeuN");
run("Set... ", "zoom=50");
setLocation(1444, 13, 999, 999);
run("Set... ", "zoom=50");
run("Brightness/Contrast...");
waitForUser("Ajust brightest. Click OK when done");
run("Threshold...");
waitForUser("Ajust threshold. Click OK when done");
run("Subtract Background...");
run("Remove Outliers...");
run("Make Binary");
run("Options...");
run("Options...");

//automatically count Cell number and CLOCK+ cell number
selectWindow("DAPI");
run("Analyze Particles...", "display add");
selectWindow("Results");
 run("Close");
nROIs = roiManager("count");
Cell = nROIs;
//s = 0;
//o = 1;
for (i=0; i<nROIs; i++) {
	selectWindow("DAPI");
	roiManager("Select", i);
	selectWindow("blackboard");
	run("Restore Selection");
	run("Add...", "value=255");
	run("Select None");
	imageCalculator("Subtract create", "blackboard","CLOCK_mask");
	run("Measure");
	selectWindow("blackboard");
	run("Multiply...", "value=0.000");
	selectWindow("Results");
	dif = getResult("Mean", 0);
	run("Close");
	if (dif > 0){
		selectWindow("DAPI_CLOCK");
		run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	    selectWindow("CLOCK+Cell");
	    run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	}
	selectWindow("Result of blackboard");
	close ();
}
selectWindow("ROI Manager");
run("Close");
selectWindow("DAPI_CLOCK");
run("Find Maxima...", "noise=10 output=Count");
CLOCK_Cell = getResult("Count", 0);
selectWindow("Results");
run("Close");

//automatically count oligos number
selectWindow("Olig2");
run("Analyze Particles...", "display add");
selectWindow("Results");
 run("Close");
nROIs = roiManager("count");
for (i=0; i<nROIs; i++) {
	selectWindow("Olig2");
	roiManager("Select", i);
	selectWindow("blackboard");
	run("Restore Selection");
	run("Add...", "value=255");
	run("Select None");
	imageCalculator("Subtract create", "blackboard","DAPI_mask");
	run("Measure");
	selectWindow("blackboard");
	run("Multiply...", "value=0.000");
	selectWindow("Results");
	dif = getResult("Mean", 0);
	run("Close");
	if (dif > 0){
		selectWindow("Olig2_DAPI");
		run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	}
	selectWindow("Result of blackboard");
	close ();
}
selectWindow("ROI Manager");
run("Close");
selectWindow("Olig2_DAPI");
run("Find Maxima...", "noise=10 output=Count");
Oligodendrocyte = getResult("Count", 0);
selectWindow("Results");
run("Close");

//automatically count CLOCK+ oligos number
selectWindow("Olig2_DAPI");
run("Analyze Particles...", "display add");
selectWindow("Results");
 run("Close");
nROIs = roiManager("count");
for (i=0; i<nROIs; i++) {
	selectWindow("Olig2_DAPI");
	roiManager("Select", i);
	selectWindow("blackboard");
	run("Restore Selection");
	run("Add...", "value=255");
	run("Select None");
	imageCalculator("Subtract create", "blackboard","CLOCK_mask");
	run("Measure");
	selectWindow("blackboard");
	run("Multiply...", "value=0.000");
	selectWindow("Results");
	dif = getResult("Mean", 0);
	run("Close");
	if (dif > 0){
		selectWindow("Olig2_DAPI_CLOCK");
		run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	    selectWindow("CLOCK+Oligo");
	    run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	}
	selectWindow("Result of blackboard");
	close ();
}
selectWindow("Olig2_DAPI_CLOCK");
run("Find Maxima...", "noise=10 output=Count");
CLOCK_Oligodendrocyte = getResult("Count", 0);
selectWindow("Results");
run("Close");
selectWindow("ROI Manager");
run("Close");

//automatically count neuron number
selectWindow("NeuN");
run("Analyze Particles...", "display add");
selectWindow("Results");
 run("Close");
nROIs = roiManager("count");
for (i=0; i<nROIs; i++) {
	selectWindow("NeuN");
	roiManager("Select", i);
	selectWindow("blackboard");
	run("Restore Selection");
	run("Add...", "value=255");
	run("Select None");
	imageCalculator("Subtract create", "blackboard","DAPI_mask");
	run("Measure");
	selectWindow("blackboard");
	run("Multiply...", "value=0.000");
	selectWindow("Results");
	dif = getResult("Mean", 0);
	run("Close");
	if (dif > 0){
		selectWindow("NeuN_DAPI");
		run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	}
	selectWindow("Result of blackboard");
	close ();
}
selectWindow("ROI Manager");
run("Close");
selectWindow("NeuN_DAPI");
run("Find Maxima...", "noise=10 output=Count");
Neuron = getResult("Count", 0);
selectWindow("Results");
run("Close");

//automatically count CLOCK+ neuron number
selectWindow("NeuN_DAPI");
run("Analyze Particles...", "display add");
selectWindow("Results");
 run("Close");
nROIs = roiManager("count");
for (i=0; i<nROIs; i++) {
	selectWindow("NeuN_DAPI");
	roiManager("Select", i);
	selectWindow("blackboard");
	run("Restore Selection");
	run("Add...", "value=255");
	run("Select None");
	imageCalculator("Subtract create", "blackboard","CLOCK_mask");
	run("Measure");
	selectWindow("blackboard");
	run("Multiply...", "value=0.000");
	selectWindow("Results");
	dif = getResult("Mean", 0);
	run("Close");
	if (dif > 0){
		selectWindow("NueN_DAPI_CLOCK");
		run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	    selectWindow("CLOCK+Neuron");
	    run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	}
	selectWindow("Result of blackboard");
	close ();
}
selectWindow("NueN_DAPI_CLOCK");
run("Find Maxima...", "noise=10 output=Count");
CLOCK_Neuron = getResult("Count", 0);
selectWindow("Results");
run("Close");

//save data and images
setResult("ID", 0, name);
setResult("Cell#", 0, Cell);
setResult("CLK_Cell#", 0, CLOCK_Cell);
setResult("Neuron#", 0, Neuron);
setResult("CLK_Neuron#", 0, CLOCK_Neuron);
setResult("Oligo#", 0, Oligodendrocyte);
setResult("CLK_Oligo#", 0, CLOCK_Oligodendrocyte);
saveAs("Results",dir+name+".csv");
run("Close");
selectWindow("CLOCK+Oligo");
saveAs("Tiff",dir+name+"_CLOCK+Oligo");
selectWindow("CLOCK+Neuron");
saveAs("Tiff",dir+name+"_CLOCK+Neuron");
selectWindow("CLOCK+Cell");
saveAs("Tiff",dir+name+"_CLOCK+Cell");
selectWindow("CLOCK");

//close windows
while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
```

-----
