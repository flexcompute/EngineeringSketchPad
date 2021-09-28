// wv.js implements functions for wv
// written by John Dannenhoffer and Bob Haimes


function postMessage(mesg) {
    // alert("in postMessage(" + mesg + ")");

    var botm = document.getElementById("botmframe");

    var text = document.createTextNode(mesg);
    botm.insertBefore(text, botm.lastChild);

    var br = document.createElement("br");
    botm.insertBefore(br, botm.lastChild);

    br.scrollIntoView();
}


function resizeFrames() {

    // get the size of the client (minus 20 to account for possible scrollbars)
    var body = document.getElementById("mainBody");
    var bodyWidth  = body.clientWidth  - 20;
    var bodyHeight = body.clientHeight - 20;

    // get the elements associated with the frames and the canvas
    var left = document.getElementById("leftframe");
    var rite = document.getElementById("riteframe");
    var botm = document.getElementById("botmframe");
    var keyf = document.getElementById("keyframe");
    var canv = document.getElementById(wv.canvasID);
    if (wv.canvasKY != undefined)
      var canf = document.getElementById(wv.canvasKY);

    // compute and set the widths of the frames
    //    (do not make leftframe larger than 250px)
    var leftWidth = Math.round(0.25 * bodyWidth);
    if (leftWidth > 250)   leftWidth = 250;
    var riteWidth = bodyWidth - leftWidth;
    var canvWidth = riteWidth - 20;
    var keyfWidth = leftWidth - 20;

    left.style.width = leftWidth + "px";
    rite.style.width = riteWidth + "px";
    botm.style.width = riteWidth + "px";
    keyf.style.width = leftWidth + "px";
    canv.style.width = canvWidth + "px";
    canv.width       = canvWidth;
    if (wv.canvasKY != undefined) {
      canf.style.width = keyfWidth + "px";
      canf.width       = keyfWidth;
    }

    // compute and set the heights of the frames
    //    (do not make botm frame larger than 200px)
    var botmHeight = Math.round(0.20 * bodyHeight);
    if (botmHeight > 200)   botmHeight = 200;
    var  topHeight = bodyHeight - botmHeight;
    var canvHeight =  topHeight - 25;
    var keyfHeight = botmHeight - 25;

    left.style.height =  topHeight + "px";
    rite.style.height =  topHeight + "px";
    botm.style.height = botmHeight + "px";
    keyf.style.height = botmHeight + "px";
    canv.style.height = canvHeight + "px";
    canv.height       = canvHeight;
    if (wv.canvasKY != undefined) {
      canf.style.height = keyfHeight + "px";
      canf.height       = keyfHeight;
      // force a key redraw
      wv.drawKey = 1;
    }
}

//
// Event Handlers
//

function getCursorXY(e)
{
    if (!e) var e = event;

    wv.cursorX  = e.clientX;
    wv.cursorY  = e.clientY;
    wv.cursorX -= wv.offLeft + 1;
    wv.cursorY  = wv.height - wv.cursorY + wv.offTop + 1;

    wv.modifier = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


function getMouseDown(e)
{
    if (!e) var e = event;

    wv.startX   = e.clientX;
    wv.startY   = e.clientY;
    wv.startX  -= wv.offLeft + 1;
    wv.startY   = wv.height - wv.startY + wv.offTop + 1;

    wv.dragging = true;
    wv.button   = e.button;

    wv.modifier = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


function getMouseUp(e)
{
    wv.dragging = false;
}

//
// callback when the mouse wheel is rolled in canvas in main mode
//
function getMouseRoll(e) {
    if (e) {

        // zoom in
        if        (e.deltaY > 0) {
            wv.mvMatrix.scale(1.1, 1.1, 1.1);
            wv.scale *= 1.1;
            wv.sceneUpd = 1;

        // zoom out
        } else if (e.deltaY < 0) {
            wv.mvMatrix.scale(0.9, 0.9, 0.9);
            wv.scale *= 0.9;
            wv.sceneUpd = 1;
        }
    }
};

function mouseLeftCanvas(e)
{
    if (wv.dragging) {
        wv.dragging = false;
    }
}


function getKeyPress(e)
{
    if (!e) var e = event;

    wv.keyPress = e.charCode;
}


//
// call-back functions
//

function wvInitUI()
{

    // set up extra storage for matrix-matrix multiplies
    wv.uiMatrix = new J3DIMatrix4();

                                  // ui cursor variables
    wv.cursorX  = -1;              // current cursor position
    wv.cursorY  = -1;
    wv.keyPress = -1;              // last key pressed
    wv.startX   = -1;              // start of dragging position
    wv.startY   = -1;
    wv.button   = -1;              // button pressed
    wv.modifier =  0;              // modifier (shift,alt,cntl) bitflag
    wv.offTop   =  0;              // offset to upper-left corner of the canvas
    wv.offLeft  =  0;
    wv.dragging = false;           // true during drag operation
    wv.picking  =  0;              // keycode of command that turns picking on
//  wv.txtInit                     // message sent at text protocol initialize
//  wv.pick                        // set to 1 to turn picking on
//  wv.picked                      // sceneGraph object that was picked
//  wv.locate                      // set to 1 to turn on locating
//  wv.located                     // coordinates (local) that were pointed at
//  wv.centerV                     // set to 1 to put wv into view centering mode
//  wv.sceneUpd                    // set to 1 when scene needs rendering
//  wv.sgUpdate                    // is 1 if the sceneGraph has been updated

    var canvas = document.getElementById(wv.canvasID);
      canvas.addEventListener('mousemove',  getCursorXY,     false);
      canvas.addEventListener('mousedown',  getMouseDown,    false);
      canvas.addEventListener('mouseup',    getMouseUp,      false);
      canvas.addEventListener('wheel',      getMouseRoll,    false);
      canvas.addEventListener('mouseout',   mouseLeftCanvas, false);
    document.addEventListener('keypress',   getKeyPress,     false);
}


function wvUpdateUI()
{

    // special code for delayed-picking mode
    if (wv.picking > 0) {

        // if something is picked, post a message
        if (wv.picked !== undefined) {

            // second part of 'q' operation
            if (wv.picking == 113) {
                postMessage("Picked: " + wv.picked.gprim);
                wv.socketUt.send("Picked: " + wv.picked.gprim);
            }

            wv.picked  = undefined;
            wv.picking = 0;
            wv.pick    = 0;

        // abort picking on mouse motion
        } else if (wv.dragging) {
            postMessage("Picking aborted");

            wv.picking = 0;
            wv.pick    = 0;
        }

        wv.keyPress = -1;
        wv.dragging = false;
    }
  
    // if something is located, post a message
    if (wv.located != undefined) {
        postMessage("Located: " + wv.located[0] + " " + wv.located[1] +
                    " " + wv.located[2]);
        wv.socketUt.send("Located: " + wv.located[0] + " " + wv.located[1] +
                         " " + wv.located[2]);
        wv.located = undefined;
        wv.locate  = 0;
    }

    // if the tree has not been created but the scene graph (possibly) exists...
    if (wv.sgUpdate == 1 && (wv.sceneGraph !== undefined)) {

        // ...count number of primitives in the scene graph
        var count = 0;
        for (var gprim in wv.sceneGraph) {

            // parse the name
            var matches = gprim.split(" ");

            var ibody = Number(matches[1]);
            if        (matches[2] == "Face") {
                var iface = matches[3];
            } else if (matches[2] == "Loop") {
                var iloop = matches[3];
            } else if (matches[2] == "Edge") {
                var iedge = matches[3];
            } else {
                alert("unknown type: " + matches[2]);
                continue;
            }

            // determine if Body does not exists
            var knode = -1;
            for (var jnode = 1; jnode < myTree.name.length; jnode++) {
                if (myTree.name[jnode] == "Body " + ibody) {
                    knode = jnode;
                }
            }

            // if Body does not exist, create it and its Face, Loop, and Edge
            //    subnodes now
            var kface, kloop, kedge;
            if (knode < 0) {
                postMessage("Processing Body " + ibody);

                myTree.addNode(0, "Body " + ibody, "*");
                knode = myTree.name.length - 1;

                myTree.addNode(knode, "__Faces", "*");
                kface = myTree.name.length - 1;

                myTree.addNode(knode, "__Loops", "*");
                kloop = myTree.name.length - 1;

                myTree.addNode(knode, "__Edges", "*");
                kedge = myTree.name.length - 1;

            // otherwise, get pointers to the face-group and loop-group nodes
            } else {
                kface = myTree.child[knode];
                kloop = kface + 1;
                kedge = kloop + 1;
            }

            // make the tree node
            if        (matches[2] == "Face") {
                myTree.addNode(kface, "____face " + iface, gprim);
            } else if (matches[2] == "Loop") {
                myTree.addNode(kloop, "____loop " + iloop, gprim);
            } else if (matches[2] == "Edge") {
                myTree.addNode(kedge, "____edge " + iedge, gprim);
            }

            count++;
        }

        // if we had any primitives, we are assuming that we have all of
        //    them, so build the tree and remember that we have
        //    built the tree
        if (count > 0) {
            myTree.build();
            wv.sgUpdate = 0;
        }
    }

    // deal with key presses
    if (wv.keyPress != -1) {

        // '?' -- help
        if (wv.keyPress ==  63) {
            postMessage("C - make tessellation coarser");
            postMessage("F - make tessellation finer");
            postMessage("N - next scalar");
            postMessage("L - change scalar limits");
            postMessage("q - query at cursor");
            postMessage("l - ~XYZ  at cursor");
            postMessage("m - set view matrix");
            postMessage("M - current view matrix");
            postMessage("x - view from -X direction");
            postMessage("X - view from +X direction");
            postMessage("y - view from -Y direction");
            postMessage("Y - view from +Y direction");
            postMessage("z - view from -Z direction");
            postMessage("Z - view from +Z direction");
            postMessage("* - center view at cursor");
            postMessage("? - get help");

        // 'C' -- make tessellation coarser
        } else if (wv.keyPress == 67) {
            postMessage("Retessellating...");
            wv.socketUt.send("coarser");

        // 'F' -- make tessellation finer
        } else if (wv.keyPress == 70) {
            postMessage("Retessellating...");
            wv.socketUt.send("finer");
          
        // 'M' -- dump view matrix
        } else if (wv.keyPress == 77) {
            var mVal = wv.mvMatrix.getAsArray();
            postMessage(" View Matrix (scale = " + wv.scale + "):");
            postMessage("  "+mVal[ 0]+" "+mVal[ 1]+" "+mVal[ 2]+" "+mVal[ 3]);
            postMessage("  "+mVal[ 4]+" "+mVal[ 5]+" "+mVal[ 6]+" "+mVal[ 7]);
            postMessage("  "+mVal[ 8]+" "+mVal[ 9]+" "+mVal[10]+" "+mVal[11]);
            postMessage("  "+mVal[12]+" "+mVal[13]+" "+mVal[14]+" "+mVal[15]);
      
        // 'N' -- next scalar
        } else if (wv.keyPress == 78) {
            wv.socketUt.send("next");
      
        // 'L' -- scalar limits
        } else if (wv.keyPress == 76) {
            wv.socketUt.send("limits");
          
        // 'l' -- locate
        } else if (wv.keyPress == 108) {
            if (wv.locate == 1) {
                wv.locate = 0;
            } else {
                wv.locate = 1;
            }
        
        // 'm' -- set view matrix
        } else if (wv.keyPress == 109) {
            var mVal  = wv.mvMatrix.getAsArray();
            wv.scale  = parseFloat(prompt("Enter scale", wv.scale));
            var strng = prompt("Matrix Pos  0- 3",
                               mVal[ 0]+" "+mVal[ 1]+" "+mVal[ 2]+" "+mVal[ 3]);
            var split = strng.split(' ');
            var len   = split.length;
            if (len > 4) len = 4;
            for (var m =  0; m <    len; m++) mVal[m] = parseFloat(split[m   ]);
            strng = prompt("Matrix Pos  4- 7",
                               mVal[ 4]+" "+mVal[ 5]+" "+mVal[ 6]+" "+mVal[ 7]);
            split = strng.split(' ');
            len   = split.length;
            if (len > 4) len = 4;
            for (var m =  4; m <  4+len; m++) mVal[m] = parseFloat(split[m- 4]);
            strng = prompt("Matrix Pos  8-11",
                               mVal[ 8]+" "+mVal[ 9]+" "+mVal[10]+" "+mVal[11]);
            split = strng.split(' ');
            len   = split.length;
            if (len > 4) len = 4;
            for (var m =  8; m <  8+len; m++) mVal[m] = parseFloat(split[m- 8]);
            strng = prompt("Matrix Pos 12-15",
                               mVal[12]+" "+mVal[13]+" "+mVal[14]+" "+mVal[15]);
            split = strng.split(' ');
            len   = split.length;
            if (len > 4) len = 4;
            for (var m = 12; m < 12+len; m++) mVal[m] = parseFloat(split[m-12]);
            wv.mvMatrix.load(mVal);
            wv.sceneUpd = 1;

        // 'q' -- query at cursor
        } else if (wv.keyPress == 113) {
            wv.picking  = 113;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'x' -- view from -X direction
        } else if (wv.keyPress == 120) {
            wv.mvMatrix.makeIdentity();
            wv.mvMatrix.rotate(+90, 0,1,0);
            wv.sceneUpd = 1;

        // 'X' -- view from +X direction
        } else if (wv.keyPress ==  88) {
            wv.mvMatrix.makeIdentity();
            wv.mvMatrix.rotate(-90, 0,1,0);
            wv.sceneUpd = 1;

        // 'y' -- view from +Y direction
        } else if (wv.keyPress == 121) {
            wv.mvMatrix.makeIdentity();
            wv.mvMatrix.rotate(-90, 1,0,0);
            wv.sceneUpd = 1;

        // 'Y' -- view from +Y direction
        } else if (wv.keyPress ==  89) {
            wv.mvMatrix.makeIdentity();
            wv.mvMatrix.rotate(+90, 1,0,0);
            wv.sceneUpd = 1;

        // 'z' -- view from +Z direction
        } else if (wv.keyPress ==  122) {
            wv.mvMatrix.makeIdentity();
            wv.mvMatrix.rotate(180, 1,0,0);
            wv.sceneUpd = 1;

        // 'Z' -- view from +Z direction
        } else if (wv.keyPress ==  90) {
            wv.mvMatrix.makeIdentity();
            wv.sceneUpd = 1;
          
        // '*' -- center view at cursor
        } else if (wv.keyPress == 42) {
            wv.centerV = 1;

        } else {
            postMessage("'" + String.fromCharCode(wv.keyPress)
                        + "' is not defined.  Use '?' for help.");
        }
    }

    wv.keyPress = -1;

    // UI is in screen coordinates (not object)
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();

    // deal with mouse movement
    if (wv.dragging) {

        // cntrl is down (rotate)
        if (wv.modifier == 4) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0;
            var angleY = -(wv.startX - wv.cursorX) / 4.0;
            if ((angleX != 0.0) || (angleY != 0.0)) {
              wv.mvMatrix.rotate(angleX, 1,0,0);
              wv.mvMatrix.rotate(angleY, 0,1,0);
              wv.sceneUpd = 1;
            }

        // alt-shift is down (rotate)
        } else if (wv.modifier == 3) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0;
            var angleY = -(wv.startX - wv.cursorX) / 4.0;
            if ((angleX != 0.0) || (angleY != 0.0)) {
              wv.mvMatrix.rotate(angleX, 1,0,0);
              wv.mvMatrix.rotate(angleY, 0,1,0);
              wv.sceneUpd = 1;
            }

        // alt is down (spin)
        } else if (wv.modifier == 2) {
            var xf = wv.startX - wv.width  / 2;
            var yf = wv.startY - wv.height / 2;

            if ((xf != 0.0) || (yf != 0.0)) {
                var theta1 = Math.atan2(yf, xf);
                xf = wv.cursorX - wv.width  / 2;
                yf = wv.cursorY - wv.height / 2;

                if ((xf != 0.0) || (yf != 0.0)) {
                    var dtheta = Math.atan2(yf, xf)-theta1;
                    if (Math.abs(dtheta) < 1.5708) {
                        var angleZ = 128*(dtheta)/3.1415926;
                        wv.mvMatrix.rotate(angleZ, 0,0,1);
                        wv.sceneUpd = 1;
                    }
                }
            }

        // shift is down (zoom)
        } else if (wv.modifier == 1) {
            if (wv.cursorY != wv.startY) {
              var scale = Math.exp((wv.cursorY - wv.startY) / 512.0);
              wv.mvMatrix.scale(scale, scale, scale);
              wv.scale   *= scale;
              wv.sceneUpd = 1;
            }

        // no modifier (translate)
        } else {
            var transX = (wv.cursorX - wv.startX) / 256.0;
            var transY = (wv.cursorY - wv.startY) / 256.0;
            if ((transX != 0.0) || (transY != 0.0)) {
              wv.mvMatrix.translate(transX, transY, 0.0);
              wv.sceneUpd = 1;
            }
        }

        wv.startX = wv.cursorX;
        wv.startY = wv.cursorY;
    }
}


function wvServerMessage(text)
{
    postMessage(" Server Message: " + text);
}


function wvServerDown()
{
    alert("The server has terminated.  You will need to restart the server"
         +" and this browser in order to continue.");
}
