// Slugs.js implements functions for the Static Legacy Unstructured Geometry System
// written by John Dannenhoffer and Bob Haimes

// Copyright (C) 2010/2020  John F. Dannenhoffer, III (Syracuse University)
//
// This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//    MA  02110-1301  USA

// functions expected by wv
//    wvInitUI()
//    wvUpdateUI()
//    wvServerMessage(text)
//    wvUpdateCanvas(gl)
//    wvServerDown()

// functions associated with button presses
//    cmdTest()
//    cmdUndo()
//    cmdHelp()
//    cmdHome()
//    cmdLeft()
//    cmdRite()
//    cmdBotm()
//    cmdTop()
//    cmdIn()
//    cmdOut()

// functions associated with mouse clicks and context-menu selections
//    automaticLinksButton(e)
//    bridgeToPointButton(e)
//    colorTrianglesButton(e)
//    flattenColorButton(e)
//    deleteTriangleButton(e)
//    fillHoleButton(e)
//    generateEgadsButton(e)
//    joinPointsButton(e)
//    linkToPointButton(e)
//    markCreasesButton(e)
//    pickPointButton(e)
//    scribeToPointButton(e)
//    writeEgadsButton(e)
//    writeStlButton(e)
//    addCommentButton(e)
//    writeStpButton(e)

// functions associated with the mouse in the canvas
//    getCursorXY(e)
//    getMouseDown(e)
//    getMouseUp(e)
//    getMouseRoll(e)
//    mouseLeftCancas(e)
//    getKeyPress(e)
//    getKeyDown(e)
//    getKeyUp(e)

// functions associated with the mouse in the key window
//    setKeyLimits(e)

// functions associated with toggling display settings
//    toggleViz(e)
//    toggleGrd(e)
//    toggleTrn(e)
//    toggleOri(e)

// functions associated with a Tree in the treefrm
//    Tree(doc, treeId) - constructor
//    TreeAddNode(iparent, name, gprim, click, cmenu, prop1, valu1, cbck1, prop2, valu2, cbck2, prop3, valu3, cbck3)
//    TreeBuild()
//    TreeClear()
//    TreeContract(inode)
//    TreeExpand(inode)
//    TreeProp(inode, iprop, onoff)
//    TreeUpdate()

// helper functions
//    postMessage(mesg)
//    resizeFrames()
//    unhighlightColumn1()

"use strict";


//
// callback when the user interface is to be initialized (called by wv-render.js)
//
function wvInitUI()
{
    // alert("wvInitUI()");

    // set up extra storage for matrix-matrix multiplies
    wv.uiMatrix   = new J3DIMatrix4();
    wv.saveMatrix = new J3DIMatrix4(wv.mvMatrix);

                                   // ui cursor variables
    wv.cursorX   = -1;             // current cursor position
    wv.cursorY   = -1;
    wv.keyPress  = -1;             // last key pressed
    wv.keyCode   = -1;
    wv.startX    = -1;             // start of dragging position
    wv.startY    = -1;
    wv.button    = -1;             // button pressed
    wv.modifier  =  0;             // modifier (shift,alt,cntl) bitflag
    wv.flying    =  1;             // flying multiplier (do not set to 0)
    wv.offTop    =  0;             // offset to upper-left corner of the canvas
    wv.offLeft   =  0;
    wv.dragging  =  false;         // true during drag operation
    wv.picking   =  0;             // keycode of command that turned picking on
    wv.locating  =  0;             // keycode of command that turned locating on
    wv.focus     = [0, 0, 0, 1];   // focus data needed in locating
    wv.initialize= 0;              //  0 need to initialize
                                   // >0 waiting for initialization (request already made)
    wv.menuEvent = undefined;      // event associated with click in Tree
    wv.server    = undefined;      // string passed back from "identify;"
    wv.loLimit   = -2;             // lower limit in key
    wv.upLimit   = +2;             // upper limit in key
    wv.curMode   =  0;             //-1 to disable buttons, etc.
                                   // 0 show WebViewer in canvas
    wv.ncolr     = 0;              // maximum color number
    wv.scale     =  1;             // scale factor for axes
//  wv.centerV                     // set to 1 to center view and rotation
//  wv.pick                        // set to 1 to turn picking on
//  wv.locate                      // set to 1 to turn locating on
//  wv.sceneGraph                  // pointer to sceneGraph
//  wv.picked                      // sceneGraph object that was picked
//  wv.located                     // sceneGraph object that was located
//  wv.sceneUpd                    // should be set to 1 to re-render scene
//  wv.sgUpdate                    // =1 if the sceneGraph has been updated
//  wv.canvasID                    // ID of main canvas
//  wv.canvasKY                    // ID of key  canvas
//  wv.drawKey                     // =1 if key window needs to be redrawn
//  wv.socketUt.send(text)         // function to send text to server
//  wv.plotAttrs                   // plot attributes

    var canvas = document.getElementById(wv.canvasID);
      canvas.addEventListener('mousemove',  getCursorXY,     false);
      canvas.addEventListener('mousedown',  getMouseDown,    false);
      canvas.addEventListener('mouseup',    getMouseUp,      false);
      canvas.addEventListener('wheel',      getMouseRoll,    false);
      canvas.addEventListener('mouseout',   mouseLeftCanvas, false);

    var keycan = document.getElementById(wv.canvasKY)
      keycan.addEventListener('mouseup',    setKeyLimits,    false);

    document.addEventListener('keypress',   getKeyPress,     false);
    document.addEventListener('keydown',    getKeyDown,      false);
    document.addEventListener('keyup',      getKeyUp,        false);
}


//
// callback when the user interface should be updated (called by wv-render.js)
//
function wvUpdateUI()
{
    // alert("wvUpdateUI()");

    // special code for delayed-picking mode
    if (wv.picking > 0) {

        // if something is picked, post a message
        if (wv.picked !== undefined) {

            // second part of '^' operation
            if (wv.picking == 94) {
                var mesg = "Picked: "+wv.picked.gprim;
                postMessage(mesg);
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

    // special code for delayed-locating mode
    if (wv.locating > 0) {

        // if something is located, post a message
        if (wv.located !== undefined) {

            // second part of '@' operation
            if (wv.locating == 64) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Located: x="+xloc.toFixed(4)+", y="+yloc.toFixed(4)+
                                   ", z="+zloc.toFixed(4));

            // second part of 'b' operation
            } else if (wv.locating ==  98) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Bridge to Point at x="+xloc.toFixed(4)+
                                             ", y="+yloc.toFixed(4)+
                                             ", z="+zloc.toFixed(4));

                wv.socketUt.send("bridgeToPoint;"+xloc.toPrecision(8)+";"+
                                                  yloc.toPrecision(8)+";"+
                                                  zloc.toPrecision(8)+";");

            // second part of 'c' operation
            } else if (wv.locating ==  99) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                var icolr = wv.ncolr+1;
                for (var jcolr = 1; jcolr <= wv.ncolr; jcolr++) {
                    var available = 1;
                    for (var gprim in wv.sceneGraph) {

                        var matches = gprim.split(" ");

                        if (matches[0] == "Color" && matches[1] == jcolr) {
                            available = 0;
                            break;
                        }
                    }

                    if (available == 1) {
                        icolr = jcolr;
                        break;
                    }
                }

                icolr = prompt('Enter color number', icolr);

                postMessage("Color Triangle at x="+xloc.toFixed(4)+
                                            ", y="+yloc.toFixed(4)+
                                            ", z="+zloc.toFixed(4));

                wv.socketUt.send("colorTriangles;"+xloc.toPrecision(8)+";"+
                                                   yloc.toPrecision(8)+";"+
                                                   zloc.toPrecision(8)+";"+icolr+";");

            // second part of 'd' operation
            } else if (wv.locating == 100) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Delete Triangle at x="+xloc.toFixed(4)+
                                             ", y="+yloc.toFixed(4)+
                                             ", z="+zloc.toFixed(4));

                wv.socketUt.send("deleteTriangle;"+xloc.toPrecision(8)+";"+yloc.toPrecision(8)+";"+zloc.toPrecision(8)+";");

            // second part of 'f' operation
            } else if (wv.locating == 102) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("fill hole at x="+xloc.toFixed(4)+
                                       ", y="+yloc.toFixed(4)+
                                       ", z="+zloc.toFixed(4));

                wv.socketUt.send("fillHole;"+xloc.toPrecision(8)+";"+
                                             yloc.toPrecision(8)+";"+
                                             zloc.toPrecision(8)+";");

            // second part of 'i' operation
            } else if (wv.locating == 105) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                wv.socketUt.send("identifyPoint;"+xloc.toPrecision(8)+";"+
                                                  yloc.toPrecision(8)+";"+
                                                  zloc.toPrecision(8)+";");

            // second part of 'F' operation
            } else if (wv.locating ==  46) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                var tol = prompt("Enter tolerance:", 0.001);

                if (tol !== null) {
                    wv.socketUt.send("flattenColor;"+xloc.toPrecision(8)+";"+
                                                     yloc.toPrecision(8)+";"+
                                                     zloc.toPrecision(8)+";"+tol+";");
                } else {
                    postMessage("Illegal tolerance");
                }

            // second part of 'j' operation
            } else if (wv.locating == 106) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Join Points at x="+xloc.toFixed(4)+
                                         ", y="+yloc.toFixed(4)+
                                         ", z="+zloc.toFixed(4));

                wv.socketUt.send("joinPoints;"+xloc.toPrecision(8)+";"+
                                               yloc.toPrecision(8)+";"
                                              +zloc.toPrecision(8)+";");

            // second part of 'k' operation
            } else if (wv.locating == 107) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Kut Triangles at x="+xloc.toFixed(4)+
                                           ", y="+yloc.toFixed(4)+
                                           ", z="+zloc.toFixed(4));

                var icolr = prompt("Enter color (or -1 for all):", "-1");
                if (icolr !== null) {
                    var type = prompt("Enter 'YZ' or 'ZX' or 'XY' or 'X' or 'Y' or 'Z' or 'data'", "YZ");
                    if (type !== null) {
                        if        (type == 'YZ' || type == 'yz' || type == 'ZY' || type == 'zy') {
                            wv.socketUt.send("cutTriangles;"+icolr+";0;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else if (type == 'XZ' || type == 'xz' || type == 'ZX' || type == 'zx') {
                            wv.socketUt.send("cutTriangles;"+icolr+";1;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else if (type == 'XY' || type == 'xy' || type == 'YX' || type == 'yx') {
                            wv.socketUt.send("cutTriangles;"+icolr+";2;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else if (type == 'X' || type == 'x') {
                            wv.socketUt.send("cutTriangles;"+icolr+";3;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else if (type == 'Y' || type == 'y') {
                            wv.socketUt.send("cutTriangles;"+icolr+";4;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else if (type == 'Z' || type == 'z') {
                            wv.socketUt.send("cutTriangles;"+icolr+";5;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else if (type == "data") {
                            wv.socketUt.send("cutTriangles;"+icolr+";6;"
                                             +xloc.toPrecision(8)+";"
                                             +yloc.toPrecision(8)+";"
                                             +zloc.toPrecision(8)+";");
                        } else {
                            alert("Illegal cut type");
                        }
                    }
                }

            // second part of 'l' operation
            } else if (wv.locating == 108) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Link Point at x="+xloc.toFixed(4)+
                                        ", y="+yloc.toFixed(4)+
                                        ", z="+zloc.toFixed(4));

                wv.socketUt.send("linkToPoint;"+xloc.toPrecision(8)+";"+
                                                yloc.toPrecision(8)+";"+
                                                zloc.toPrecision(8)+";");

            // second part of 'p' operation
            } else if (wv.locating == 112) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Pick Point at x="+xloc.toFixed(4)+
                                        ", y="+yloc.toFixed(4)+
                                        ", z="+zloc.toFixed(4));

                wv.socketUt.send("pickPoint;"+xloc.toPrecision(8)+";"+
                                              yloc.toPrecision(8)+";"+
                                              zloc.toPrecision(8)+";");

            // second part of 's' operation
            } else if (wv.locating == 115) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Scribe to x="+xloc.toFixed(4)+
                                    ", y="+yloc.toFixed(4)+
                                    ", z="+zloc.toFixed(4));

                wv.socketUt.send("scribeToPoint;"+xloc.toPrecision(8)+";"+
                                                  yloc.toPrecision(8)+";"+
                                                  zloc.toPrecision(8)+";");

            // second part of 't' operation
            } else if (wv.locating == 116) {
                var xloc = wv.focus[3] * wv.located[0] + wv.focus[0];
                var yloc = wv.focus[3] * wv.located[1] + wv.focus[1];
                var zloc = wv.focus[3] * wv.located[2] + wv.focus[2];

                postMessage("Identify Triangle  at x="+xloc.toFixed(4)+
                                                ", y="+yloc.toFixed(4)+
                                                ", z="+zloc.toFixed(4));

                wv.socketUt.send("identifyTriangle;"+xloc.toPrecision(8)+";"+
                                                     yloc.toPrecision(8)+";"+
                                                     zloc.toPrecision(8)+";");
            }

            wv.located  = undefined;
            wv.locating = 0;
            wv.locate   = 0;

        // abort locating on mouse motion
        } else if (wv.dragging) {
            postMessage("Locating aborted");

            wv.locating = 0;
            wv.locate   = 0;
        }

        wv.keyPress = -1;
        wv.dragging = false;
    }

    /* if the server has not been identified, ask it for its identity */
    if (wv.initialize > 0) {
        wv.initialize--;
    } else if (wv.initialize == 0) {
        try {
            wv.socketUt.send("identify;");
            wv.initialize = 6000;
        } catch (x) {
            // could not send command
        }
    }

    // if the scene graph has been updated, (re-)build the Tree
    if (wv.sgUpdate == 1) {

        if (wv.sceneGraph === undefined) {
            alert("wv.sceneGraph is undefined --- but we need it");
        }

        // clear previous Nodes from the Tree
        myTree.clear();

        // put the group headers into the Tree
        myTree.addNode(0, "Commands", "", null, "");                   // node 1
        myTree.addNode(0, "Display", "",  null, "");                   // node 2

        myTree.addNode(2, "\u00a0\u00a0Triangles", "", null, "",
                       "Viz", "on",  toggleViz,
                       "Grd", "on",  toggleGrd,
                       "Trn", "off", toggleTrn);                       // node 3
        myTree.addNode(1, "\u00a0\u00a0Cleanup",       "", null, "");  // node 4
        myTree.addNode(1, "\u00a0\u00a0Coloring",      "", null, "");  // node 5
        myTree.addNode(1, "\u00a0\u00a0GenerateEgads",  "", generateEgadsButton, "*", "(E)");
        myTree.addNode(1, "\u00a0\u00a0Miscellaneous", "", null, "");  // node 7


        // put the Commands into the Tree
        myTree.addNode(4, "\u00a0\u00a0\u00a0\u00a0Pick Point",        "", pickPointButton,      "*", "(p)");
        myTree.addNode(4, "\u00a0\u00a0\u00a0\u00a0Join Points",       "", joinPointsButton,     "*", "(j)");
        myTree.addNode(4, "\u00a0\u00a0\u00a0\u00a0Bridge to Point",   "", bridgeToPointButton,  "*", "(b)");
        myTree.addNode(4, "\u00a0\u00a0\u00a0\u00a0Fill hole",         "", fillHoleButton,       "*", "(f)");
        myTree.addNode(4, "\u00a0\u00a0\u00a0\u00a0Delete Triangle",   "", deleteTriangleButton, "*", "(d)");

        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Pick Point",        "", pickPointButton,      "*", "(p)");
        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Link to Point",     "", linkToPointButton,    "*", "(l)");
        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Scribe to Point",   "", scribeToPointButton,  "*", "(s)");
        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Mark Creases",      "", markCreasesButton,    "*", "(m)");
        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Color Triangles",   "", colorTrianglesButton, "*", "(c)");
        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Flatten Color",     "", flattenColorButton,   "*", "(F)");
        myTree.addNode(5, "\u00a0\u00a0\u00a0\u00a0Auto. Links",       "", automaticLinksButton, "*", "(a)");

        myTree.addNode(7, "\u00a0\u00a0\u00a0\u00a0Write .stl file",   "", writeStlButton,       "*"       );
        myTree.addNode(7, "\u00a0\u00a0\u00a0\u00a0Add comment",       "", addCommentButton,     "*", "(#)");

        // open the Display (by default)
        myTree.opened[1] = 1;           // commands
        myTree.opened[2] = 1;           // display
        myTree.opened[3] = 1;           // triangles
        myTree.opened[4] = 1;           // cleanup
        myTree.opened[5] = 1;           // coloring
        myTree.opened[7] = 1;           // miscellaneous

        // put the Display elements into the Tree
        for (var gprim in wv.sceneGraph) {

            // parse the name
            var matches = gprim.split(" ");

            if (matches[0] == "Color") {
                var icolr = Number(matches[1]);

                myTree.addNode(3, "\u00a0\u00a0\u00a0\u00a0Color "+icolr, gprim, null, "",
                               "Viz", "on",  toggleViz,
                               "Grd", "on",  toggleGrd,
                               "Trn", "off", toggleTrn);

                if (icolr > wv.ncolr) {
                    wv.ncolr = icolr;
                }

            } else if (matches[0] == "Hangs") {
                myTree.addNode(2, "\u00a0\u00a0Hangs", gprim, null, "",
                               "Viz", "on", toggleViz);

//$$$                myTree.opened[4] = 1;   // cleanup
//$$$                myTree.opened[5] = 0;   // coloring
            } else if (matches[0] == "Links") {
                myTree.addNode(2, "\u00a0\u00a0Links", gprim, null, "",
                               "Viz", "on", toggleViz);

            } else if (matches[0] == "CurPt") {
                myTree.addNode(2, "\u00a0\u00a0CurPt", gprim, null, "",
                               "Viz", "on", toggleViz);

            }
        }

        // mark that we have (re-)built the Tree
        wv.sgUpdate = 0;

        // convert the abstract Tree Nodes into an HTML table
        myTree.build();
    }

    // deal with key presses
    if (wv.keyPress != -1) {

        var myKeyPress = String.fromCharCode(wv.keyPress);

        // '?' -- help
        if (myKeyPress == "?") {
            postMessage("........................... Viewer Cursor options ...........................\n" +
                        "ctrl-h <Home> - initial view             ctrl-f        - front view          \n" +
                        "ctrl-l        - leftside view            ctrl-r        - riteside view       \n" +
                        "ctrl-t        - top view                 ctrl-b        - bottom view         \n" +
                        "ctrl-i <PgUp> - zoom in                  ctrl-o <PgDn> - zoom out            \n" +
                        "<Left>        - rotate or xlate left     <Rite>        - rotate or xlate rite\n" +
                        "<Up>          - rotate or xlate up       <Down>        - rotate or xlate down\n" +
                        ">             - save view                <             - recall view         \n" +
                        "^             - query object at cursor   @             - get coords. @ cursor\n" +
                        "!             - toggle flying mode       *             - center view @ cursor\n" +
                        ".............................. ? - get help .................................");

        // 'a' -- automatic Links
        } else if (myKeyPress == 'a') {
            automaticLinksButton(0);
    
        // 'b' -- bridge to Point
        } else if (myKeyPress == 'b') {
            wv.locating =  98;
            wv.locate   = 1;

        // 'c' -- color Triangle
        } else if (myKeyPress == 'c') {
            wv.locating =  99;
            wv.locate   = 1;

        // 'd' -- delete Triangle
        } else if (myKeyPress == 'd') {
            wv.locating = 100;
            wv.locate   = 1;

        // 'E' -- generate Egads
        } else if (myKeyPress == 'E') {
            generateEgadsButton(0);

        // 'f' -- fill hole
        } else if (myKeyPress == 'f') {
            wv.locating = 102;
            wv.locate   = 1;

        // 'F' -- flatten color
        } else if (myKeyPress == 'F') {
            wv.locating =  46;
            wv.locate   = 1;

        // 'i' -- identify Point
        } else if (myKeyPress == 'i') {
            wv.locating = 105;
            wv.locate   = 1;

        // 'j' -- join Points
        } else if (myKeyPress == 'j') {
            wv.locating = 106;
            wv.locate   = 1;

        // 'k' -- cut Triangles
        } else if (myKeyPress == 'k') {
            wv.locating = 107;
            wv.locate   = 1;

        // 'l' -- Link to Point
        } else if (myKeyPress == 'l') {
            wv.locating = 108;
            wv.locate   = 1;

        // 'm' -- mark creases
        } else if (myKeyPress == 'm') {
            markCreasesButton();

        // 'p' -- pick Point
        } else if (myKeyPress == 'p') {
            wv.locating = 112;
            wv.locate   = 1;

        // 's' -- scribe to Point
        } else if (myKeyPress == 's') {
            wv.locating = 115;
            wv.locate   = 1;

        // 't' -- identify Triangle
        } else if (myKeyPress == 't') {
            wv.locating = 116;
            wv.locate   = 1;

        // '#' -- add comment
        } else if (myKeyPress == '#') {
            addCommentButton(0);

        // '^' -- query at cursor
        } else if (myKeyPress == '^') {
            wv.picking  = 94;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '@' -- locate at cursor
        } else if (myKeyPress == '@') {
            wv.locating = 64;
            wv.locate   = 1;

        // '*' -- center view
        } else if (myKeyPress == '*') {
            wv.centerV = 1;
            postMessage("View being centered");

        // '!' -- toggle flying mode
        } else if (myKeyPress == "!") {
            if (wv.flying <= 1) {
                postMessage("Turning flying mode ON");
                wv.flying = 10;
            } else {
                postMessage("Turning flying mode OFF");
                wv.flying = 1;
            }

        // '>' -- save view
        } else if (myKeyPress == ">") {
            postMessage("Saving current view");
            wv.saveMatrix.load(wv.mvMatrix);
            wv.sceneUpd = 1;

        // '<' -- recall view
        } else if (myKeyPress == "<") {
            postMessage("Restoring saved view");
            wv.mvMatrix.load(wv.saveMatrix);
            wv.sceneUpd = 1;

        // '<esc>' -- not used
//        } else if (wv.keyPress == 0 && wv.keyCode == 27) {
//            postMessage("<Esc> is not supported.  Use '?' for help");

        // '<Home>' -- initial view
        } else if (wv.keyPress == 0 && wv.keyCode == 36) {
            wv.mvMatrix.makeIdentity();
            wv.scale    = 1;
            wv.sceneUpd = 1;

        // '<end>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 35) {
            postMessage("<End> is not supported.  Use '?' for help");

        // '<PgUp>' -- zoom in
        } else if (wv.keyPress == 0 && wv.keyCode == 33) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(2.0, 2.0, 2.0);
                wv.scale *= 2.0;
            } else {
                wv.mvMatrix.scale(1.25, 1.25, 1.25);
                wv.scale *= 1.25;
            }
            wv.sceneUpd = 1;

        // '<PgDn>' -- zoom out
        } else if (wv.keyPress == 0 && wv.keyCode == 34) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(0.5, 0.5, 0.5);
                wv.scale *= 0.5;
            } else {
                wv.mvMatrix.scale(0.8, 0.8, 0.8);
                wv.scale *= 0.8;
            }
            wv.sceneUpd = 1;

        // '<Delete>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 46) {
            postMessage("<Delete> is not supported.  Use '?' for help");

        // '<Left>' -- rotate or translate object left
        } else if (wv.keyPress == 0 && wv.keyCode == 37) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(-30, 0,1,0);
                } else {
                    wv.mvMatrix.rotate( -5, 0,1,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(-0.5, 0.0, 0.0);
                } else {
                    wv.mvMatrix.translate(-0.1, 0.0, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Right>' -- rotate or translate object right
        } else if (wv.keyPress == 0 && wv.keyCode == 39) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(+30, 0,1,0);
                } else {
                    wv.mvMatrix.rotate( +5, 0,1,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(+0.5, 0.0, 0.0);
                } else {
                    wv.mvMatrix.translate(+0.1, 0.0, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Up>' -- rotate or translate object up
        } else if (wv.keyPress == 0 && wv.keyCode == 38) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(-30, 1,0,0);
                } else {
                    wv.mvMatrix.rotate( -5, 1,0,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(0.0, +0.5, 0.0);
                } else {
                    wv.mvMatrix.translate(0.0, +0.1, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // '<Down>' -- rotate or translate object down
        } else if (wv.keyPress == 0 && wv.keyCode == 40) {
            if (wv.flying == 1) {
                if (wv.modifier == 0) {
                    wv.mvMatrix.rotate(+30, 1,0,0);
                } else {
                    wv.mvMatrix.rotate( +5, 1,0,0);
                }
            } else {
                if (wv.modifier == 0) {
                    wv.mvMatrix.translate(0.0, -0.5, 0.0);
                } else {
                    wv.mvMatrix.translate(0.0, -0.1, 0.0);
                }
            }
            wv.sceneUpd = 1;

        // 'ctrl-h' - initial view (same as <Home>)
        } else if ((wv.keyPress == 104 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   8 && wv.modifier == 4 && wv.keyCode ==  8)   ) {
//$$$            wv.mvMatrix.makeIdentity();
//$$$            wv.scale    = 1;
//$$$            wv.sceneUpd = 1;
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-i' - zoom in (same as <PgUp> without shift)
        } else if ((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) {
//$$$            wv.mvMatrix.scale(2.0, 2.0, 2.0);
//$$$            wv.scale   *= 2.0;
//$$$            wv.sceneUpd = 1;
            cmdIn();
            wv.keyPress = -1;
            return;

        // '+' - zoom in (same as <PgUp> without shift)
        } else if (wv.keyPress ==  43 && wv.modifier == 1) {
            cmdIn();
            wv.keyPress = -1;
            return;

        // 'ctrl-o' - zoom out (same as <PgDn> without shift)
        } else if ((wv.keyPress == 111 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  15 && wv.modifier == 4 && wv.keyCode == 15)   ) {
//$$$            wv.mvMatrix.scale(0.5, 0.5, 0.5);
//$$$            wv.scale   *= 0.5;
//$$$            wv.sceneUpd = 1;
            cmdOut();
            wv.keyPress = -1;
            return;

        // '-' - zoom out (same as <PgDn> without shift)
        } else if (wv.keyPress ==  45 && wv.modifier == 0) {
            cmdOut();
            wv.keyPress = -1;
            return;

        // 'ctrl-f' - front view (same as <Home>)
        } else if ((wv.keyPress == 102 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   6 && wv.modifier == 4 && wv.keyCode ==  6)   ) {
//$$$            wv.mvMatrix.makeIdentity();
//$$$            wv.scale    = 1;
//$$$            wv.sceneUpd = 1;
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-r' - riteside view
        } else if ((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) {
//$$$            wv.mvMatrix.makeIdentity();
//$$$            wv.mvMatrix.rotate(-90, 0,1,0);
//$$$            wv.sacle    = 1;
//$$$            wv.sceneUpd = 1;
            cmdRite();
            wv.keyPress = -1;
            return;

        // 'ctrl-l' - leftside view
        } else if ((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) {
//$$$            wv.mvMatrix.makeIdentity();
//$$$            wv.mvMatrix.rotate(+90, 0,1,0);
//$$$            wv.scale    = 1;
//$$$            wv.sceneUpd = 1;
            cmdLeft();
            wv.keyPress = -1;
            return;

        // 'ctrl-t' - top view
        } else if ((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) {
//$$$            wv.mvMatrix.makeIdentity();
//$$$            wv.mvMatrix.rotate(+90, 1,0,0);
//$$$            wv.scale    = 1;
//$$$            wv.sceneUpd = 1;
            cmdTop();
            wv.keyPress = -1;
            return;

        // 'ctrl-b' - bottom view
        } else if ((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) {
//$$$            wv.mvMatrix.makeIdentity();
//$$$            wv.mvMatrix.rotate(-90, 1,0,0);
//$$$            wv.scale    = 1;
//$$$            wv.sceneUpd = 1;
            cmdBotm();
            wv.keyPress = -1;
            return;

        // NOP
        } else if (wv.keyPress == 0 && wv.modifier == 0) {

//$$$        // unknown command
//$$$        } else {
//$$$            postMessage("Unknown command (keyPress="+wv.keyPress
//$$$                        +", modifier="+wv.modifier
//$$$                        +", keyCode="+wv.keyCode+").  Use '?' for help");
        }

        wv.keyPress = -1;
    }

    // UI is in screen coordinates (not object)
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();

    // deal with mouse movement
    if (wv.dragging) {

        // cntrl is down (rotate)
        if (wv.modifier == 4) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0 / wv.flying;
            var angleY = -(wv.startX - wv.cursorX) / 4.0 / wv.flying;
            if ((angleX != 0.0) || (angleY != 0.0)) {
                wv.mvMatrix.rotate(angleX, 1,0,0);
                wv.mvMatrix.rotate(angleY, 0,1,0);
                wv.sceneUpd = 1;
            }

        // alt-shift is down (rotate)
        } else if (wv.modifier == 3) {
            var angleX =  (wv.startY - wv.cursorY) / 4.0 / wv.flying;
            var angleY = -(wv.startX - wv.cursorX) / 4.0 / wv.flying;
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
                    var dtheta = Math.atan2(yf, xf) - theta1;
                    if (Math.abs(dtheta) < 1.5708) {
                        var angleZ = 128*(dtheta) / 3.1415926 / wv.flying;
                        wv.mvMatrix.rotate(angleZ, 0,0,1);
                        wv.sceneUpd = 1;
                    }
                }
            }

        // shift is down (zoom)
        } else if (wv.modifier == 1) {
            if (wv.cursorY != wv.startY) {
                var scale = Math.exp((wv.cursorY - wv.startY) / 512.0 / wv.flying);
                wv.mvMatrix.scale(scale, scale, scale);
                wv.scale   *= scale;
                wv.sceneUpd = 1;
            }

        // no modifier (translate)
        } else {
            var transX = (wv.cursorX - wv.startX) / 256.0 / wv.flying;
            var transY = (wv.cursorY - wv.startY) / 256.0 / wv.flying;
            if ((transX != 0.0) || (transY != 0.0)) {
                wv.mvMatrix.translate(transX, transY, 0.0);
                wv.sceneUpd = 1;
            }
        }

        // if not flying, then update the start coordinates
        if (wv.flying <= 1) {
            wv.startX = wv.cursorX;
            wv.startY = wv.cursorY;
        }
    }
}


//
// callback when a (text) message is received from the server (called by wv-socket.js)
//
function wvServerMessage(text)
{
    // alert("wvServerMessage(text="+text+")");

    // if it is a null message, do nothing
    if (text.length <= 1) {

    // if it starts with "identify;" post a message */
    } else if (text.substring(0,9) == "identify;") {
        if (wv.server === undefined) {
            wv.server = text.substring(9,text.length-2);
            postMessage("Slugs has been initialized and is attached to server '"+wv.server+"'");
        }

    // if it starts with "sgFocus;" store the scene graph focus data
    } else if (text.substring(0,8) == "sgFocus;") {
        if (text.length > 10) {
            wv.focus = JSON.parse(text.substring(8,text.length-1));
        } else {
            wv.focus = [];
        }

    // if it contains "ERROR::" post the error
    } else if (text.indexOf("ERROR::") >= 0) {
        postMessage(text.substring(text.indexOf("ERROR::")));
        alert("Error encountered (see messages)");

    // default is to post the message
    } else {
        postMessage(text);
    }
}


//
// callback when the server goes down either on abort or normal closure (called by wv-socket.js)
//
function wvServerDown()
{
    // deactivate the buttons
    wv.curMode = -1;

    postMessage("The server has terminated or network connection has been lost.\n"
          +"Slugs will be working off-line.");

    // turn the background of the message window pink
    var botm = document.getElementById("brframe");
    botm.style.backgroundColor = "#FFAFAF";
}


//
// callback used to put axes on the canvas
//
function wvUpdateCanvas(gl)
{
    // Construct the identity as projection matrix and pass it in
    wv.mvpMatrix.load(wv.perspectiveMatrix);
    wv.mvpMatrix.setUniform(gl, wv.u_modelViewProjMatrixLoc, false);

    var mv   = new J3DIMatrix4();
    var mVal = wv.mvMatrix.getAsArray();

    mVal[ 3] = 0.0;
    mVal[ 7] = 0.0;
    mVal[11] = 0.0;
    mv.load(mVal);
    mv.scale(1.0/wv.scale, 1.0/wv.scale, 1.0/wv.scale);
    mv.invert();
    mv.transpose();

    // define location of axes in space
    var x    = -1.5 * wv.width / wv.height;
    var y    = -1.5;
    var z    =  0.9;
    var pdel =  0.15;     // length of axes in positive direction
    var mdel =  0.0;      // length of axes in negative direction

    // set up coordinates for axes
    mVal = mv.getAsArray();

    var vertices = new Float32Array(18);
    vertices[ 0] = x - mdel*mVal[ 0];   vertices[ 1] = y - mdel*mVal[ 1];   vertices[ 2] = z - mdel*mVal[ 2];
    vertices[ 3] = x + pdel*mVal[ 0];   vertices[ 4] = y + pdel*mVal[ 1];   vertices[ 5] = z + pdel*mVal[ 2];
    vertices[ 6] = x - mdel*mVal[ 4];   vertices[ 7] = y - mdel*mVal[ 5];   vertices[ 8] = z - mdel*mVal[ 6];
    vertices[ 9] = x + pdel*mVal[ 4];   vertices[10] = y + pdel*mVal[ 5];   vertices[11] = z + pdel*mVal[ 6];
    vertices[12] = x - mdel*mVal[ 8];   vertices[13] = y - mdel*mVal[ 9];   vertices[14] = z - mdel*mVal[10];
    vertices[15] = x + pdel*mVal[ 8];   vertices[16] = y + pdel*mVal[ 9];   vertices[17] = z + pdel*mVal[10];

    // set up colors for the axes
    var colors = new Uint8Array(18);
    colors[ 0] = 255;   colors[ 1] =   0;   colors[ 2] =   0;
    colors[ 3] = 255;   colors[ 4] =   0;   colors[ 5] =   0;
    colors[ 6] =   0;   colors[ 7] = 255;   colors[ 8] =   0;
    colors[ 9] =   0;   colors[10] = 255;   colors[11] =   0;
    colors[12] =   0;   colors[13] =   0;   colors[14] = 255;
    colors[15] =   0;   colors[16] =   0;   colors[17] = 255;

    // draw the axes
    gl.lineWidth(3);
    gl.disableVertexAttribArray(2);
    gl.uniform1f(wv.u_wLightLoc, 0.0);

    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(0);

    var cbuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
    gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
    gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.enableVertexAttribArray(1);

    gl.drawArrays(gl.LINES, 0, 6);
    gl.deleteBuffer(buffer);
    gl.deleteBuffer(cbuf);
    gl.uniform1f(wv.u_wLightLoc, 1.0);
}


//
// callback when "testButton" is pressed (called by ESP.html)
//
function cmdTest () {
    alert("in cmdTest()");

}


//
// callback when "undoButton" is pressed (called by ESP.html)
//
function cmdUndo () {
//    alert("in cmdUndo()");

    wv.socketUt.send("undo;");
}


//
// callback when "helpButton" is pressed (called by ESP.html)
//
function cmdHelp () {
//    alert("cmdHelp()");

    // open help in another tab
    window.open("Slugs-help.html");
}


//
// callback when "homeButton" is pressed (called by ESP.html)
//
function cmdHome()
{
    wv.mvMatrix.makeIdentity();
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
}


//
// callback when "leftButton" is pressed (called by ESP.html)
//
function cmdLeft()
{
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(+90, 0,1,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
}


//
// callback when "riteButton" is pressed (called by ESP.html)
//
function cmdRite()
{
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(-90, 0,1,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
}


//
// callback when "botmButton" is pressed (called by ESP.html)
//
function cmdBotm()
{
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(-90, 1,0,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
}


//
// callback when "topButton" is pressed (called by ESP.html)
//
function cmdTop()
{
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(+90, 1,0,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
}


//
// callback when "inButton" is pressed (called by ESP.html)
//
function cmdIn()
{
    wv.mvMatrix.load(wv.uiMatrix);
    wv.mvMatrix.scale(2.0, 2.0, 2.0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale   *= 2.0;
    wv.sceneUpd = 1;
}


//
// callback when "outButton" is pressed (called by ESP.html)
//
function cmdOut()
{
    wv.mvMatrix.load(wv.uiMatrix);
    wv.mvMatrix.scale(0.5, 0.5, 0.5);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale   *= 0.5;
    wv.sceneUpd = 1;
}


//
// callback when "AutomaticLinks" is pressed in Tree
//
function automaticLinksButton(e) {
//    alert("in automaticLinksButton(e"+e+")");

    wv.socketUt.send("automaticLinks;");
}


//
// callback when "BridgeToPoint" is pressed in Tree
//
function bridgeToPointButton(e) {
//    alert("in bridgeToPointButton(e="+e+")");

    alert("> Set CurrentPoint with \"p\" command.\n"+
          "> Point to a Triangle with a Hanging side and press \"b\" command.\n"+
          "-----\n"+
          "< New Triangle formed between them.\n"+
          "< CurrentPoint is unset.\n"+
          "-----\n"+
          "The \"b\" command can be undone");
}


//
// callback when "ColorTriangles" is pressed in Tree
//
function colorTrianglesButton(e) {
//    alert("in colorTrianglesButton(e="+e+")");

    alert("> Point to a Triangle and press \"c\" command.\n"+
          "> Enter color index.\n"+
          "-----\n"+
          "< All contigous Triangles (not across a Link) are colored\n"+
          "< CurrentPoint is unset.\n"+
          "-----\n"+
          "The \"c\" command can be undone.");
}


//
// callback when "FlattenColor" is pressed in Tree
//
function flattenColorButton(e) {
//    alert("in flattenColorButton(e="+e+")");

    alert("> Point to a Triangle and press \"F\" command.\n"+
          "-----\n"+
          "The \"F\" command can be undone.");
}


//
// callback when "DeleteTriangle" is pressed in Tree
//
function deleteTriangleButton(e) {
//    alert("in deleteTriangleButton(e="+e+")");

    alert("> Point to a Triangle and press \"d\" command.\n"+
          "-----\n"+
          "< Selected Triangle is deleted\n"+
          "< Hanging sides are updated\n"+
          "< CurrentPoint is unset.\n"+
          "-----\n"+
          "The \"d\" command can be undone.");
}


//
// callback when "FillHole" is pressed in Tree
//
function fillHoleButton(e) {
//    alert("in fillHoleButton(e="+e+")");

    alert("> Point to Triangle with at least one Hanging side and press \"f\" command.\n"+
          "-----\n"+
          "< Triangles are created to fill the hole.\n"+
          "< CurrentPoint is unset.\n"+
          "-----\n"+
          "The \"f\" command can be undone.");
}


//
// callback when "GenerateEgads" is pressed in Tree
//
function generateEgadsButton(e) {
//    alert("in generateEgadsButton(e="+e+")");

    var filename = prompt("Enter filename of EGADS file:");
    if (filename !== null) {
        var len = filename.length;
        if (len < 6) {
            filename += ".egads";
        } else if (filename.substring(len-6) != ".egads") {
            filename += ".egads";
        }

        var ncp = prompt("Enter number of control points:", 7);
        if (ncp !== null) {
            wv.socketUt.send("generateEgads;"+filename+";"+ncp+";");
        }
    }
}


//
// callback when "JoinPoints" is pressed in Tree
//
function joinPointsButton(e) {
//    alert("in joinPointsButton(e="+e+")");

    alert("> Set the CurrentPoint by using the \"p\" command.\n"+
          "> Pick another Point to join and press \"j\" command.\n"+
          "-----\n"+
          "< New Point is created at average of these Points.\n"+
          "< CurrentPoint is set to new Point.\n"+
          "-----\n"+
          "The \"j\" command can be undone");
}


//
// callback when "LinkToPoint" is pressed in Tree
//
function linkToPointButton(e) {
//    alert("in linkToPointButton(e="+e+")");

    alert("> Select the CurrentPoint by using the \"p\" command\n."+
          "> Pick new Point and press the (letter) \"l\" command.\n"+
          "-----\n"+
          "< Links are created on shortest path between the two Points.\n"+
          "< CurrentPoint is moved to Point at \"l\" command.\n"+
          "-----\n"+
          "The \"l\" command can be undone.");
}


//
// callback when "MarkCreases" is pressed in Tree
//
function markCreasesButton(e) {
//    alert("in markCreasesButton(e="+e+")");

    var angdeg = prompt("Enter minimum crease angle (deg):", 45);
    if (angdeg !==null) {
        wv.socketUt.send("markCreases;"+angdeg+";");
    }
}


//
// callback when "PickPoint" is pressed in Tree
//
function pickPointButton(e) {
//    alert("in pickPointButton(e="+e+")");

    alert("> Select a Point and press the \"p\" command.\n"+
          "-----\n"+
          "< CurrentPoint is set.\n"+
          "-----\n"+
          "The \"p\" command does not affect the ability to undo a previous command.");
}


//
// callback when "ScribeToPoint" is pressed in Tree
//
function scribeToPointButton(e) {
//    alert("in scribePointToButton(e="+e+")");

    alert("> Select the CurrentPoint by using the \"p\" command\n."+
          "> Pick new Point and press \"s\" command.\n"+
          "-----\n"+
          "< Links are created on shortest path between the two Points.\n"+
          "< CurrentPoint is moved to Point at \"s\" command.\n"+
          "-----\n"+
          "The \"s\" command can be undone.");
}


//
// callback when "WriteEgads" is pressed in Tree
//
function writeEgadsButton(e) {
//    alert("in writeEgadsButton(e="+e+")");

    alert("Write .egads file is not implemented.");
}


//
// callback when "WriteStl" is pressed in Tree
//
function writeStlButton(e) {
//    alert("in writeStlButton(e="+e+")");

    // add .stl if not present
    var filename = prompt("Enter .stl filename:");
    if (filename.indexOf(".stl") < 0) {
        filename = filename.concat(".stl");
    }

    postMessage("Write .stl file \""+filename+"\"");

    wv.socketUt.send("writeStlFile;"+filename+";");
}


//
// callback when "AddComment" is pressed in Tree
//
function addCommentButton(e) {
//    alert("in addCommentButton(e="+e+")");

    // add .stl if not present
    var comment = prompt("Enter comment:");

    postMessage("Add comment \""+comment+"\"");

    wv.socketUt.send("addComment;"+comment+";");
}


//
// callback when "WriteStp" is pressed in Tree
//
function writeStpButton(e) {
//    alert("in writeStpButton(e="+e+")");

    alert("Write .stp file is not implemented.");
}


//
// callback to get the cursor location
//
function getCursorXY(e)
{
    if (!e) var e = event;

    wv.cursorX  =  e.clientX - wv.offLeft             - 1;
    wv.cursorY  = -e.clientY + wv.offTop  + wv.height + 1;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


//
// callback when any mouse is pressed in canvas
//
function getMouseDown(e)
{
    if (!e) var e = event;

    wv.startX   =  e.clientX - wv.offLeft             - 1;
    wv.startY   = -e.clientY + wv.offTop  + wv.height + 1;

    wv.dragging = true;
    wv.button   = e.button;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
}


//
// callback when the mouse is released
//
function getMouseUp(e)
{
    wv.dragging = false;
}


//
// callback when the mouse wheel is rolled in canvas
//
function getMouseRoll(e)
{
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
}


//
// callback when the mouse leaves the canvas
//
function mouseLeftCanvas(e)
{
    if (wv.dragging) {
        wv.dragging = false;
    }
}


//
// callback when a key is pressed
//
function getKeyPress(e)
{
    if (wv.curMode != 0) {
        return;
    }

//    if (!e) var e = event;

    wv.keyPress = e.charCode;
    wv.keyCode  = e.keyCode;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;

    // make sure that Firefox does not stop websockets if <esc> is pressed
    if (e.charCode == 0) {
        e.preventDefault();
    }
}


//
// callback when an arrow... or shift is pressed (needed for Chrome)
//
function getKeyDown(e) {
    if (e.charCode == 0) {
        if (e.keyCode == 33 ||          // PgUp
            e.keyCode == 34 ||          // PgDn
            e.keyCode == 36 ||          // Home
            e.keyCode == 37 ||          // Left
            e.keyCode == 38 ||          // Up
            e.keyCode == 39 ||          // Right
            e.keyCode == 40   ) {       // Down
            wv.keyCode  = e.keyCode;
            wv.keyPress = 0;
        } else if (e.keyCode == 16) {   // Shift
            wv.modifier = 1;
        }
    }
}


//
// callback when a shift is released (needed for Chrome)
//
function getKeyUp(e) {
    if (e.charCode == 0 && e.keyCode == 16) {
        wv.modifier = 0;
    }
}


//
// callback when the mouse is pressed in key window
//
function setKeyLimits(e)
{

    // get new limits
    var templo = prompt("Enter new lower limit", wv.loLimit);
    if (isNaN(templo)) {
        alert("Lower limit must be a number");
        return;
    }

    var tempup = prompt("Enter new upper limit", wv.upLimit);
    if (isNaN(tempup)) {
        alert("Upper limit must be a number");
        return;
    }

    if (templo != wv.loLimit || tempup != wv.uplimit) {
        wv.loLimit = templo;
        wv.upLimit = tempup;

        // send the limits back to the server
        wv.socketUt.send("setLims;"+wv.loLimit+";"+wv.upLimit+";");
    }
}


//
// callback to toggle Viz property
//
function toggleViz(e) {
    // alert("in toggleViz(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first.");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Viz property
    if        (myTree.valu1[inode] == "off") {
        myTree.prop(inode, 1, "on");
    } else if (myTree.valu1[inode] == "on") {
        myTree.prop(inode, 1, "off");
    } else {
        alert("illegal Viz property:"+myTree.valu1[inode]);
        return;
    }
}


//
// callback to toggle Grd property
//
function toggleGrd(e) {
    // alert("in toggleGrd(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first.");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Grd property
    if        (myTree.valu2[inode] == "off") {
        myTree.prop(inode, 2, "on");
    } else if (myTree.valu2[inode] == "on") {
        myTree.prop(inode, 2, "off");
    } else {
        alert("illegal Grd property:"+myTree.valu2[inode]);
        return;
    }
}


//
// callback to toggle Trn property
//
function toggleTrn(e) {
    // alert("in toggleTrn(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first.");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Trn property
    if        (myTree.valu3[inode] == "off") {
        myTree.prop(inode, 3, "on");
    } else if (myTree.valu3[inode] == "on") {
        myTree.prop(inode, 3, "off");
    } else {
        alert("illegal Trn property:"+myTree.valu3[inode]);
        return;
    }
}


//
// callback to toggle Ori property
//
function toggleOri(e) {
    // alert("in toggleOri(e="+e+")");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first.");
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Ori property
    if        (myTree.valu3[inode] == "off") {
        myTree.prop(inode, 3, "on");
    } else if (myTree.valu3[inode] == "on") {
        myTree.prop(inode, 3, "off");
    } else {
        alert("illegal Ori property:"+myTree.valu3[inode]);
        return;
    }
}


//
// constructor for a Tree
//
function Tree(doc, treeId) {
    // alert("in Tree(doc="+doc+", treeId="+treeId+")");

    // remember the document
    this.document = doc;
    this.treeId   = treeId;

    // arrays to hold the Nodes
    this.name   = new Array();
    this.gprim  = new Array();
    this.click  = new Array();
    this.cmenu  = new Array();
    this.parent = new Array();
    this.child  = new Array();
    this.next   = new Array();
    this.nprop  = new Array();
    this.opened = new Array();

    this.prop1  = new Array();
    this.valu1  = new Array();
    this.cbck1  = new Array();
    this.prop2  = new Array();
    this.valu2  = new Array();
    this.cbck2  = new Array();
    this.prop3  = new Array();
    this.valu3  = new Array();
    this.cbck3  = new Array();

    // initialize Node=0 (the root)
    this.name[  0] = "**root**";
    this.gprim[ 0] = "";
    this.click[ 0] = null;
    this.cmenu[ 0] = "";
    this.parent[0] = -1;
    this.child[ 0] = -1;
    this.next[  0] = -1;
    this.nprop[ 0] =  0;
    this.prop1[ 0] = "";
    this.valu1[ 0] = "";
    this.cbck1[ 0] = null;
    this.prop2[ 0] = "";
    this.valu2[ 0] = "";
    this.cbck2[ 0] = null;
    this.prop3[ 0] = "";
    this.valu3[ 0] = "";
    this.cbck3[ 0] = null;
    this.opened[0] = +1;

    // add methods
    this.addNode  = TreeAddNode;
    this.expand   = TreeExpand;
    this.contract = TreeContract;
    this.prop     = TreeProp;
    this.clear    = TreeClear;
    this.build    = TreeBuild;
    this.update   = TreeUpdate;
}


//
// add a Node to the Tree
//
function TreeAddNode(iparent, name, gprim, click, cmenu, prop1, valu1, cbck1,
                     prop2, valu2, cbck2, prop3, valu3, cbck3)
{
    // alert("in TreeAddNode(iparent="+iparent+", name="+name+", gprim="+gprim+", click="+click+
    //                       ", cmenu="+cmenu+", prop1="+prop1+", valu1="+valu1+",cbck1="+cbck1+
    //                       ", prop2="+prop2+", valu2="+valu2+", cbck2="+cbck2+", prop3="+prop3+
    //                       ", valu3="+valu3+", cbck3="+cbck3+")");

    // validate the input
    if (iparent < 0 || iparent >= this.name.length) {
        alert("iparent="+iparent+" is out of range");
        return;
    }

    // find the next Node index
    var inode = this.name.length;

    // store this Node's values
    this.name[  inode] = name;
    this.gprim[ inode] = gprim;
    this.click[ inode] = click;
    this.cmenu[ inode] = cmenu;
    this.parent[inode] = iparent;
    this.child[ inode] = -1;
    this.next[  inode] = -1;
    this.nprop[ inode] =  0;
    this.opened[inode] =  0;

    // store the properties
    if (prop1 !== undefined) {
        this.nprop[inode] = 1;
        this.prop1[inode] = prop1;
        this.valu1[inode] = valu1;
        this.cbck1[inode] = cbck1;
    }

    if (prop2 !== undefined) {
        this.nprop[inode] = 2;
        this.prop2[inode] = prop2;
        this.valu2[inode] = valu2;
        this.cbck2[inode] = cbck2;
    }

    if (prop3 !== undefined) {
        this.nprop[inode] = 3;
        this.prop3[inode] = prop3;
        this.valu3[inode] = valu3;
        this.cbck3[inode] = cbck3;
    }

    // if the parent does not have a child, link this
    //    new Node to the parent
    if (this.child[iparent] < 0) {
        this.child[iparent] = inode;

    // otherwise link this Node to the last parent's child
    } else {
        var jnode = this.child[iparent];
        while (this.next[jnode] >= 0) {
            jnode = this.next[jnode];
        }

        this.next[jnode] = inode;
    }
}


//
// build the Tree (ie, create the html table from the Nodes)
//
function TreeBuild() {
    // alert("in TreeBuild()");

    var doc = this.document;

    // if the table already exists, delete it and all its children (3 levels)
    var thisTable = doc.getElementById(this.treeId);
    if (thisTable) {
        var child1 = thisTable.lastChild;
        while (child1) {
            var child2 = child1.lastChild;
            while (child2) {
                var child3 = child2.lastChild;
                while (child3) {
                    child2.removeChild(child3);
                    child3 = child2.lastChild;
                }
                child1.removeChild(child2);
                child2 = child1.lastChild;
            }
            thisTable.removeChild(child1);
            child1 = thisTable.lastChild;
        }
        thisTable.parentNode.removeChild(thisTable);
    }

    // build the new table
    var newTable = doc.createElement("table");
    newTable.setAttribute("id", this.treeId);
    doc.getElementById("treefrm").appendChild(newTable);

    // traverse the Nodes using depth-first search
    var inode = 1;
    while (inode > 0) {

        // table row "node"+inode
        var newTR = doc.createElement("TR");
        newTR.setAttribute("id", "node"+inode);
        newTable.appendChild(newTR);

        // table data "node"+inode+"col1"
        var newTDcol1 = doc.createElement("TD");
        newTDcol1.setAttribute("id", "node"+inode+"col1");
        newTDcol1.className = "fakelinkon";
        newTR.appendChild(newTDcol1);

        var newTexta = doc.createTextNode("");
        newTDcol1.appendChild(newTexta);

        // table data "node"+inode+"col2"
        var newTDcol2 = doc.createElement("TD");
        newTDcol2.setAttribute("id", "node"+inode+"col2");
        if (this.cmenu[inode] != "") {
            newTDcol2.className = "fakelinkcmenu";
        }
        newTR.appendChild(newTDcol2);

        var newTextb = doc.createTextNode(this.name[inode]);
        newTDcol2.appendChild(newTextb);

        var name = this.name[inode].replace(/\u00a0/g, "");

        var ibrch = 0;
        for (var jbrch = 0; jbrch < brch.length; jbrch++) {
            if (brch[jbrch].name == name) {
                if (brch[jbrch].ileft == -2) {
                    newTDcol2.className = "errorTD";
                }
                break;
            }
        }

        // table data "node"+inode+"col3"
        if (this.nprop[inode] > 0) {
            var newTDcol3 = doc.createElement("TD");
            newTDcol3.setAttribute("id", "node"+inode+"col3");
            if (this.cbck1[inode] != "") {
                newTDcol3.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol3);

            if (this.nprop[inode] == 1) {
                newTDcol3.setAttribute("colspan", "3");
            }

            var newTextc = doc.createTextNode(this.prop1[inode]);
            newTDcol3.appendChild(newTextc);
        }

        // table data "node:+inode+"col4"
        if (this.nprop[inode] > 1) {
            var newTDcol4 = doc.createElement("TD");
            newTDcol4.setAttribute("id", "node"+inode+"col4");
            if (this.cbck2[inode] != "") {
                newTDcol4.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol4);

            if (this.nprop[inode] == 2) {
                newTDcol4.setAttribute("colspan", "2");
            }

            var newTextd = doc.createTextNode(this.prop2[inode]);
            newTDcol4.appendChild(newTextd);
        }

        // table data "node:+inode+"col5"
        if (this.nprop[inode] > 2) {
            var newTDcol5 = doc.createElement("TD");
            newTDcol5.setAttribute("id", "node"+inode+"col5");
            if (this.cbck3[inode] != "") {
                newTDcol5.className = "fakelinkon";
            }
            newTR.appendChild(newTDcol5);

            var newTextd = doc.createTextNode(this.prop3[inode]);
            newTDcol5.appendChild(newTextd);
        }

        // go to next row
        if        (this.child[inode] >= 0) {
            inode = this.child[inode];
        } else if (this.next[inode] >= 0) {
            inode = this.next[inode];
        } else {
            while (inode > 0) {
                inode = this.parent[inode];
                if (this.parent[inode] == 0) {
                    newTR = doc.createElement("TR");
                    newTR.setAttribute("height", "10px");
                    newTable.appendChild(newTR);
                }
                if (this.next[inode] >= 0) {
                    inode = this.next[inode];
                    break;
                }
            }
        }
    }

    this.update();
}


//
// clear the Tree
//
function TreeClear() {
    // alert("in TreeClear()");

    // remove all but the first Node
    this.name.splice(  1);
    this.gprim.splice( 1);
    this.click.splice( 1);
    this.cmenu.splice( 1);
    this.parent.splice(1);
    this.child.splice( 1);
    this.next.splice(  1);
    this.nprop.splice( 1);
    this.opened.splice(1);

    this.prop1.splice(1);
    this.valu1.splice(1);
    this.cbck1.splice(1);
    this.prop2.splice(1);
    this.valu2.splice(1);
    this.cbck2.splice(1);
    this.prop3.splice(1);
    this.valu3.splice(1);
    this.cbck3.splice(1);

    // reset the root Node
    this.parent[0] = -1;
    this.child[ 0] = -1;
    this.next[  0] = -1;
}


//
// expand a Node in the Tree
//
function TreeContract(inode) {
    // alert("in TreeContract(inode="+inode+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    }

    // contract inode
    this.opened[inode] = 0;

    // contract all descendents of inode
    for (var jnode = 1; jnode < this.parent.length; jnode++) {
        var iparent = this.parent[jnode];
        while (iparent > 0) {
            if (iparent == inode) {
                this.opened[jnode] = 0;
                break;
            }

            iparent = this.parent[iparent];
        }
    }

    // update the display
    this.update();
}


//
// expand a Node in the Tree
//
function TreeExpand(inode) {
    // alert("in TreeExpand(inode="+inode+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    }

    // expand inode
    this.opened[inode] = 1;

    // update the display
    this.update();
}


//
// change a property of a Node
//
function TreeProp(inode, iprop, onoff) {
    // alert("in TreeProp(inode="+inode+", iprop="+iprop+", onoff="+onoff+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    } else if (onoff != "on" && onoff != "off") {
        alert("onoff="+onoff+" is not 'on' or 'off'");
        return;
    }

    // set the property for inode
    if (iprop == 1) {
        this.valu1[inode] = onoff;
    } else if (iprop == 2) {
        this.valu2[inode] = onoff;
    } else if (iprop == 3) {
        this.valu3[inode] = onoff;
    } else {
        alert("iprop="+iprop+" is not 1, 2, or 3");
        return;
    }

    // set property of all descendents of inode
    for (var jnode = 1; jnode < this.parent.length; jnode++) {
        var iparent = this.parent[jnode];
        while (iparent > 0) {
            if (iparent == inode) {
                if        (iprop == 1) {
                    this.valu1[jnode] = onoff;
                } else if (iprop == 2) {
                    this.valu2[jnode] = onoff;
                } else if (iprop == 3) {
                    this.valu3[jnode] = onoff;
                }
                break;
            }

            iparent = this.parent[iparent];
        }
    }

    this.update();
}


//
// update the Tree (after build/expension/contraction/property-set)
//
function TreeUpdate() {
//    alert("in TreeUpdate()");

    var doc = this.document;

    // traverse the Nodes using depth-first search
    for (var inode = 1; inode < this.opened.length; inode++) {
        var element = doc.getElementById("node"+inode);

        // unhide the row
        element.style.display = null;

        // hide the row if one of its parents has .opened=0
        var jnode = this.parent[inode];
        while (jnode != 0) {
            if (this.opened[jnode] == 0) {
                element.style.display = "none";
                break;
            }

            jnode = this.parent[jnode];
        }

        // if the current Node has children, set up appropriate event handler to expand/collapse
        if (this.child[inode] > 0) {
            if (this.opened[inode] == 0) {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.firstChild.nodeValue = "+";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.expand(thisNode);
                };

            } else {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.firstChild.nodeValue = "-";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.contract(thisNode);
                };
            }
        }

        if (this.click[inode] !== null) {
            var myElem = doc.getElementById("node"+inode+"col2");
            myElem.onclick = this.click[inode];
        }

        // set the class of the properties
        if (this.nprop[inode] >= 1) {
            var myElem = doc.getElementById("node"+inode+"col3");
            myElem.onclick = this.cbck1[inode];

            if (this.prop1[inode] == "Viz") {
                if (this.valu1[inode] == "off") {
                    myElem.setAttribute("class", "fakelinkoff");
                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ON;
                        wv.sceneUpd = 1;
                    }
                } else {
                    myElem.setAttribute("class", "fakelinkon");
                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ON;
                        wv.sceneUpd = 1;
                    }
                }
            }
        }

        if (this.nprop[inode] >= 2) {
            var myElem = doc.getElementById("node"+inode+"col4");
            myElem.onclick = this.cbck2[inode];

            if (this.prop2[inode] == "Grd") {
                if (this.valu2[ inode] == "off") {
                    myElem.setAttribute("class", "fakelinkoff");

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.LINES;
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.POINTS;
                        wv.sceneUpd = 1;
                    }
                } else {
                    myElem.setAttribute("class", "fakelinkon");

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.LINES;
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.POINTS;
                        wv.sceneUpd = 1;
                    }
                }
            }
        }

        if (this.nprop[inode] >= 3) {
            var myElem = doc.getElementById("node"+inode+"col5");
            myElem.onclick = this.cbck3[inode];

            if (this.prop3[inode] == "Trn") {
                if (this.valu3[ inode] == "off") {
                    myElem.setAttribute("class", "fakelinkoff");

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        wv.sceneUpd = 1;
                    }
                } else {
                    myElem.setAttribute("class", "fakelinkon");

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.TRANSPARENT;
                        wv.sceneUpd = 1;
                    }
                }
            } else if (this.prop3[inode] == "Ori") {
                if (this.valu3[ inode] == "off") {
                    myElem.setAttribute("class", "fakelinkoff");

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ORIENTATION;
                        wv.sceneUpd = 1;
                    }
                } else {
                    myElem.setAttribute("class", "fakelinkon");

                    if (this.gprim[inode] != "") {
                        wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ORIENTATION;
                        wv.sceneUpd = 1;
                    }
                }
            }
        }
    }
}


//
// post a message into the brframe
//
function postMessage(mesg) {
    // alert("in postMessage(mesg="+mesg+")");

    var botm = document.getElementById("brframe");

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);
    botm.insertBefore(pre, botm.lastChild);

    pre.scrollIntoView();
}


//
// callback when "onresize" event occurs (called by ESP.html)
//
// resize the frames (with special handling to width of tlframe and height of brframe
//
function resizeFrames() {
    // alert("resizeFrames()");

    var scrollSize;
    if (BrowserDetect.browser == "Chrome") {
        scrollSize = 24;
    } else {
        scrollSize = 20;
    }

    // get the size of the client (minus amount to account for scrollbars)
    var body = document.getElementById("mainBody");
    var bodyWidth  = body.clientWidth  - scrollSize;
    var bodyHeight = body.clientHeight - scrollSize;

    // get the elements associated with the frames and the canvas
    var topleft = document.getElementById("tlframe");
    var butnfrm = document.getElementById("butnfrm");
    var treefrm = document.getElementById("treefrm");
    var toprite = document.getElementById("trframe");
    var botleft = document.getElementById("blframe")
    var botrite = document.getElementById("brframe");
    var canv    = document.getElementById(wv.canvasID);

    // compute and set the widths of the frames
    //    (do not make tlframe larger than 250px)
    var leftWidth = Math.round(0.25 * bodyWidth);
    if (leftWidth > 250)   leftWidth = 250;
    var riteWidth = bodyWidth - leftWidth;
    var canvWidth = riteWidth - scrollSize;

    topleft.style.width = leftWidth+"px";
    butnfrm.style.width = leftWidth+"px";
    treefrm.style.width = leftWidth+"px";
    toprite.style.width = riteWidth+"px";
    botleft.style.width = leftWidth+"px";
    botrite.style.width = riteWidth+"px";
    canv.style.width    = canvWidth+"px";
    canv.width          = canvWidth;

    // compute and set the heights of the frames
    //    (do not make brframe larger than 200px)
    var botmHeight = Math.round(0.20 * bodyHeight);
    if (botmHeight > 200)   botmHeight = 200;
    var  topHeight = bodyHeight - botmHeight;
    var canvHeight =  topHeight - scrollSize - 5;
    var keyHeight  = botmHeight - 25;

    topleft.style.height =  topHeight+"px";
    treefrm.style.height = (topHeight-80)+"px";
    toprite.style.height =  topHeight+"px";
    botleft.style.height = botmHeight+"px";
    botrite.style.height = botmHeight+"px";
    canv.style.height    = canvHeight+"px";
    canv.height          = canvHeight;

    // set up canvas associated with key
    if (wv.canvasKY != undefined) {
        var canf = document.getElementById(wv.canvasKY);

        var keyfWidth     = leftWidth - 20;
        canf.style.width  = keyfWidth + "px";
        canf.width        = keyfWidth;

        var keyfHeight    = botmHeight - 25;
        canf.style.height = keyfHeight + "px";
        canf.height       = keyfHeight;

        // force a key redraw
        wv.drawKey = 1;
    }
}


//
// unhighlight background colors in column 1 of myTree
//
function unhighlightColumn1() {
    // alert("in unhighlightColumn1()");

    for (var jnode = 0; jnode < myTree.name.length; jnode++) {
        var myElem = document.getElementById("node"+jnode+"col1");
        if (myElem) {
            if (myElem.className == "currentTD" || myElem.className == "parentTD" || myElem.className == "childTD") {
                myElem.className = undefined;
            }
        }
    }
}
