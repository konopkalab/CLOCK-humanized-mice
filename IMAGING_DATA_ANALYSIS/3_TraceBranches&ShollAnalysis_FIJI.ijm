//all numbers in the script need to be modified base on the images and computer
//open image, adjust size and set image location on screen, and get the name of image
open()
name=getTitle;
name=substring(name,0,lengthOf(name)-4);
setLocation(-1900, 13, 4000, 4000);
run("Set... ", "zoom=13");
setTool("hand");

//mannually adjust brightness
run("Brightness/Contrast...");
waitForUser("Ajust brightness. Click OK when done");

//Start Simple Neurite Tracer and mannually trace branches
run("8-bit");
run("Simple Neurite Tracer", "look_for_tubeness look_for_previously_traced");
waitForUser("Trace all dendrites. Click OK when done");

//save traced data for banches
run("IJ Robot", "order=Left_Click x_point=53 y_point=73 delay=40");
run("IJ Robot", "order=Left_Click x_point=53 y_point=137 delay=40");
waitForUser("SAVE. Click OK when done");
run("IJ Robot", "order=Left_Click x_point=20 y_point=73 delay=40");
run("IJ Robot", "order=Left_Click x_point=20 y_point=144 delay=40");
waitForUser("SAVE. Click OK when done");
run("IJ Robot", "order=Left_Click x_point=53 y_point=73 delay=40");
run("IJ Robot", "order=Left_Click x_point=53 y_point=118 delay=400");

//mannually draw a line from the center of soma to the farthest tip of branch
setTool("line");
selectWindow("Paths rendered in a Stack");
run("Options...", "iterations=5 count=1 black do=Dilate");
waitForUser("Draw a radius line. Click OK when done");

//run Sholl analysis and save results
selectWindow("Paths rendered in a Stack");
run("Sholl Analysis...", "starting=5 radius_step=1 #_samples=1 integration=Mean enclosing=1 #_primary=1 infer fit linear polynomial=[Best fitting degree] most normalizer=Area directory=[]");
dir = "C:\\Users\\aaronfans007\\Desktop\\Neurite\\UPPER\\";
selectWindow("Paths rendered in a Stack");
saveAs("Jpeg",dir+name+"_traced");
selectWindow("Paths rendered in a Stack_Sholl-Profiles");
saveAs("Results", dir+name+"_sholl.csv");

//close windows
while (nImages>0) { 
        selectImage(nImages); 
        close(); 
        } 
run("IJ Robot", "order=Left_Click x_point=850 y_point=154 delay=40");
run("IJ Robot", "order=Left_Click x_point=880 y_point=154 delay=40");
run("IJ Robot", "order=Left_Click x_point=812 y_point=44 delay=40");
run("IJ Robot", "order=Left_Click x_point=974 y_point=44 delay=40");
run("IJ Robot", "order=Left_Click x_point=313 y_point=44 delay=40");
run("IJ Robot", "order=Left_Click x_point=672 y_point=282 delay=40");
