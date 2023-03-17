# LAMINATION AND LAYER THICKNESS
- [1_Lamination&LayerThickness_FIJI.ijm](1_Lamination&LayerThickness_FIJI.ijm)

```
//open image, split channel, adjust window size and location
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
```

-----
