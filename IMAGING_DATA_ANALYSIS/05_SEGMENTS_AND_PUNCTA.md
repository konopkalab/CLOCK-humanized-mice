# SEGMENTS AND PUNCTA

```
//get the list of files that will be processed in a batch and open files one by one
fileList = getFileList("C:/Users/aaronfans007/Desktop/temp/a/");
for (n = 0; n < fileList.length; n++)
{
open("C:/Users/aaronfans007/Desktop/temp/a/"+fileList[n]);
dir = "C:/Users/aaronfans007/Desktop/temp/b/";
name=getTitle;
name_f=substring(name,37,lengthOf(name)-4);

//split channels
run("Duplicate...", "title=dup duplicate");
run("Duplicate...", "title=DAPI duplicate channels=1");
selectWindow("dup");
run("Duplicate...", "title=TUJ1 duplicate channels=3");
selectWindow("dup");
run("Duplicate...", "title=vGLUT2 duplicate channels=4");

//create a whiteboard
selectWindow("DAPI");
run("Duplicate...", "title=whiteboard");
run("Make Binary");
selectWindow("whiteboard");
run("Multiply...", "value=0.000");
run("Duplicate...", "title=DAPI_final duplicate");

//calculate number of nuclei (DAPI signal) that are neurons
selectWindow("DAPI");
run("8-bit");
setAutoThreshold("Yen");
setOption("BlackBackground", false);
run("Convert to Mask");
run("Invert");
run("Analyze Particles...", "size=23-Infinity circularity=0.13-1.00 display add");
selectWindow("Results");
run("Close");
nROIs = roiManager("count");
  for (i=0; i<nROIs; i++) {  
    selectWindow("DAPI");
	roiManager("Select", i);
	selectWindow("whiteboard");
	run("Restore Selection");
	run("Add...", "value=255");
	run("Select None");
	run("Measure");
	mean_before = getResult("Mean", 0);
	selectWindow("Results");
    run("Close");
	imageCalculator("Subtract create", "whiteboard","TUJ1_mask");
	run("Measure");
	mean_after = getResult("Mean", 0);
	selectWindow("Results");
    run("Close");
    selectWindow("whiteboard");
    run("Multiply...", "value=0.000");
	if (mean_after > mean_before*0.73){
		selectWindow("DAPI_final");
		run("Restore Selection");
	    run("Add...", "value=255");
	    run("Select None");
	}
	selectWindow("Result of whiteboard");
	close ();
  }	
selectWindow("ROI Manager");
run("Close");
selectWindow("DAPI_final");
run("Analyze Particles...", "size=23-Infinity circularity=0.13-1.00 display add");
saveAs("Results",dir+name_f+"_DAPI.csv");
selectWindow("Results");
run("Close");
selectWindow("ROI Manager");
run("Close");

//apply Analyze Skeleton function to quantify segments of TUJ1 signal for the complexity of dendritic arborization
selectWindow("TUJ1");
run("Duplicate...", "title=TUJ1_mask");
selectWindow("TUJ1_mask");
run("8-bit");
setAutoThreshold("Default");
run("Convert to Mask");
selectWindow("TUJ1");
run("Subtract Background...", "rolling=50");
run("8-bit");
setAutoThreshold("Yen");
run("Convert to Mask");
run("Invert");
run("Analyze Skeleton (2D/3D)", "prune=none");
saveAs("Results",dir+name_f+"_branch.csv");
selectWindow("Results");
run("Close");

//quantify density of vGLUT2 puncta on neuronal branches (TUJ1 signal)
selectWindow("vGLUT2");
run("Subtract Background...", "rolling=50");
setAutoThreshold("Otsu");
run("Convert to Mask");
run("Invert");
selectWindow("dup");
run("Duplicate...", "title=TUJ1_mask1 duplicate channels=3");
run("Subtract Background...", "rolling=50");
run("8-bit");
setAutoThreshold("Yen");
run("Convert to Mask");
imageCalculator("Subtract create", "vGLUT2","TUJ1_mask1");
run("Find Maxima...", "noise=10 output=Count light");
saveAs("Results",dir+name_f+"_puncta.csv");
selectWindow("Results");
run("Close");

//close all windows
while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
}
```
