//open image and split channels
open()
name=getTitle;
run("Split Channels");
selectWindow("C1-"+name);
run("Blue");
run("Set... ", "zoom=10");
setLocation(0, 0, 700, 770);
run("Set... ", "zoom=10");
selectWindow("C2-"+name);
run("Green");
run("Set... ", "zoom=10");
setLocation(450, 0, 700, 770);
run("Set... ", "zoom=10");
selectWindow("C3-"+name);
run("Red");
run("Set... ", "zoom=10");
setLocation(897, 0, 700, 770);
run("Set... ", "zoom=10");
AllenSection = getNumber("Please enter a value for its corresponding Allen Atlas section",0);
selectWindow("C2-"+name);
run("Set... ", "zoom=50");
setLocation(0, 0, 700, 770);
run("Set... ", "zoom=50");
setTool("hand");
waitForUser("Please go to ROI_ScaleBar. Click OK when done");

//convert the unit from pxiel to micron
x=1;
while (x==1) {
	setTool("line");
	waitForUser("Please draw line on scale bar. Click OK when done");
	if (selectionType() == 5) {
		distance = getNumber("Please enter a value for the real distance of scale bar",1000);
		run("Set Scale...", "known=distance pixel=1 unit=micron");
		x=0;
	}
}

//measure the area of hemisphere in coronal section
selectWindow("C1-"+name);
run("Set... ", "zoom=12");
setLocation(0, 0, 700, 770);
run("Set... ", "zoom=12");
run("Brightness/Contrast...");
showMessage("Adjust brightness if neccessary");
setTool("polygon");
x=1;
while (x==1) {
	waitForUser("Please draw polygon to measure area of half section. Click OK when done");
	if (selectionType() == 2) {
		run("Measure");
		Overlay.addSelection("yellow", 7); 
		dir = getDirectory("choose a directory");
		saveAs("Tiff",dir+name+"_DAPI_area");
		x=0;
	}
}

//measure cortical thickness
selectWindow("C3-"+name);
run("Set... ", "zoom=23");
setLocation(0, 0, 1300, 769);
run("Set... ", "zoom=23");
setTool("hand");
waitForUser("Please zoom in and go to ROI. Click OK when done");
setTool("line");
for (m=1; m<4; m++) {
	selectWindow("C3-"+name);
	waitForUser("Please draw line "+m+" to measure Layer I-VI. Click OK when done");
	if (selectionType() == 5) {
		run("Measure");
		Overlay.addSelection("yellow", 13); 
		selectWindow("C2-"+name);
		run("Restore Selection");
		Overlay.addSelection("yellow", 13);
}
}

//measure layer 6
selectWindow("C3-"+name);
run("Set... ", "zoom=50");
setLocation(0, 0, 1300, 770);
run("Set... ", "zoom=50");
setTool("hand");
waitForUser("Please zoom in and go to ROI. Click OK when done");
setTool("line");
for (n=1; n<4; n++) {
	waitForUser("Please draw line "+n+" to measure Layer VI. Click OK when done");
	if (selectionType() == 5) {
		run("Measure");
		Overlay.addSelection("blue", 7); 
}
}
saveAs("Tiff",dir+name+"_FOXP2_LayerVI");

//measure upper layers
selectWindow("C2-"+name);
run("Set... ", "zoom=40");
setLocation(0, 0, 1300, 770);
run("Set... ", "zoom=40");
setTool("hand");
waitForUser("Please zoom in and go to ROI. Click OK when done");
setTool("line");
for (i=1; i<4; i++) {
	waitForUser("Please draw line "+i+" to measure Layer I-IV. Click OK when done");
	if (selectionType() == 5) {
		run("Measure");
		Overlay.addSelection("red", 13); 
}
}
for (j=1; j<4; j++) {
	waitForUser("Please draw line "+j+" to measure Layer II/IV. Click OK when done");
	if (selectionType() == 5) {
		run("Measure");
		Overlay.addSelection("blue", 7); 
}
}
saveAs("Tiff",dir+name+"_CUX1_LayerI-IV");
selectWindow(name+"_DAPI_area.tif");
close();
selectWindow(name+"_FOXP2_LayerVI.tif");
close();
selectWindow(name+"_CUX1_LayerI-IV.tif");
close();

//save data to speadsheet
area = getResult("Area",0);
Cortex1 = getResult("Length",1);
Cortex2 = getResult("Length",2);
Cortex3 = getResult("Length",3);
Ave_Cortex = (Cortex1 + Cortex2 + Cortex3)/3;
LayerVI1 = getResult("Length",4);
LayerVI2 = getResult("Length",5);
LayerVI3 = getResult("Length",6);
Ave_LayerVI = (LayerVI1 + LayerVI2 + LayerVI3)/3;
LayerUpper1 = getResult("Length",7);
LayerUpper2 = getResult("Length",8);
LayerUpper3 = getResult("Length",9);
Ave_LayerUpper = (LayerUpper1 + LayerUpper2 + LayerUpper3)/3;
LayerIItoIV1 = getResult("Length",10);
LayerIItoIV2 = getResult("Length",11);
LayerIItoIV3 = getResult("Length",12);
Ave_LayerIItoIV = (LayerIItoIV1 + LayerIItoIV2 + LayerIItoIV3)/3;
LayerI1 = LayerUpper1 - LayerIItoIV1;
LayerI2 = LayerUpper2 - LayerIItoIV2;
LayerI3 = LayerUpper3 - LayerIItoIV3;
Ave_LayerI = (LayerI1 + LayerI2 + LayerI3)/3;
R_Cortex1 = Cortex1/area;
R_Cortex2 = Cortex2/area;
R_Cortex3 = Cortex3/area;
Ave_R_Cortex = (R_Cortex1 + R_Cortex2 + R_Cortex3)/3;
R_LayerI1 = LayerI1/Cortex1;
R_LayerI2 = LayerI2/Cortex2;
R_LayerI3 = LayerI3/Cortex3;
Ave_R_LayerI = (LayerI1 + LayerI2 + LayerI3)/3;
R_LayerIItoIV1 = LayerIItoIV1/Cortex1;
R_LayerIItoIV2 = LayerIItoIV2/Cortex2;
R_LayerIItoIV3 = LayerIItoIV3/Cortex3;
Ave_R_LayerIItoIV = (R_LayerIItoIV1 + R_LayerIItoIV2 + R_LayerIItoIV3)/3;
R_LayerVI1 = LayerVI1/Cortex1;
R_LayerVI2 = LayerVI2/Cortex3;
R_LayerVI3 = LayerVI3/Cortex3;
Ave_R_LayerVI = (R_LayerVI1 + R_LayerVI2 + R_LayerVI3)/3;
selectWindow("Results");
run("Close");
setResult("ID", 0, name);
setResult("AllenSection", 0, AllenSection);
setResult("area", 0, area);
setResult("Cortex1", 0, Cortex1);
setResult("Cortex2", 0, Cortex2);
setResult("Cortex3", 0, Cortex3);
setResult("Ave_Cortex", 0, Ave_Cortex);
setResult("R_Cortex1", 0, R_Cortex1);
setResult("R_Cortex2", 0, R_Cortex2);
setResult("R_Cortex3", 0, R_Cortex3);
setResult("Ave_R_Cortex", 0, Ave_R_Cortex);
setResult("LayerVI1", 0, LayerVI1);
setResult("LayerVI2", 0, LayerVI2);
setResult("LayerVI3", 0, LayerVI3);
setResult("Ave_LayerVI", 0, Ave_LayerVI);
setResult("R_LayerVI1", 0, R_LayerVI1);
setResult("R_LayerVI2", 0, R_LayerVI2);
setResult("R_LayerVI3", 0, R_LayerVI3);
setResult("Ave_R_LayerVI", 0, Ave_R_LayerVI);
setResult("LayerIItoIV1", 0, LayerIItoIV1);
setResult("LayerIItoIV2", 0, LayerIItoIV2);
setResult("LayerIItoIV3", 0, LayerIItoIV3);
setResult("Ave_LayerIItoIV", 0, Ave_LayerIItoIV);
setResult("R_LayerIItoIV1", 0, R_LayerIItoIV1);
setResult("R_LayerIItoIV2", 0, R_LayerIItoIV2);
setResult("R_LayerIItoIV3", 0, R_LayerIItoIV3);
setResult("Ave_R_LayerIItoIV", 0, Ave_R_LayerIItoIV);
saveAs("Results",dir+name+".csv");
selectWindow("Results");
run("Close");
