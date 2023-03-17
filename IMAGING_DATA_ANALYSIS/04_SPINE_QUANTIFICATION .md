# SPINE QUANTIFICATION

```
//open, relocate, and adjust brightness of branch segment image
open()
name=getTitle;
name_f=substring(name,0,lengthOf(name)-4);
setLocation(-1900, 13, 4000, 4000);
run("Set... ", "zoom=300");
run("Brightness/Contrast...");
waitForUser("brightness/contrast. Click OK when done");

//measure the length of branch segment
setTool("polyline");
waitForUser("draw a line. Click OK when done");
run("Measure");
saveAs("Results", "C:/Users/aaronfans007/Desktop/Neurite/UPPER/branch/"+name_f+".txt");

//run SpineJ and save spine quantification results
selectWindow(name);
run("Duplicate...", "title=dup");
setLocation(-1900, 13, 4000, 4000);
run("Set... ", "zoom=300");
setTool("hand");
selectWindow(name);
run("SpineJ");
waitForUser("Click OK when done");
name=getTitle;
name_f=substring(name,0,lengthOf(name)-4);
saveAs("Text","C:/Users/aaronfans007/Desktop/Neurite/UPPER/branch/"+name_f+".xls");

//close all windows
while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
selectWindow("Log");
run("Close");
selectWindow("Results");
run("Close");
selectWindow(name_f+".xls");
run("Close");
```