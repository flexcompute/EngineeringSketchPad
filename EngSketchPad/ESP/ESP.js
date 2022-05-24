// ESP.js implements functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer and Bob Haimes

// Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
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

// interface functions that a tool can provide (provided below for tool=main)
//    launch()
//
//    cmdUndo()
//    cmdSolve()
//    cmdSave()
//    cmdQuit()
//
//    cmdHome()            button presses
//    cmdLeft()
//    cmdRite()
//    cmdBotm()
//    cmdTop()
//    cmdIn()
//    cmdOut()
//
//    mouseDown(e)         mouse events
//    mouseMove(e)
//    mouseUp(e)
//    mouseWheel(e)
//    mouseLeftCanvas(e)
//
//    keyPress(e)          key presses
//    keyDown(e)
//    keyUp(e)
//    keyPressPart1(myKeyPress)
//    keyPressPart2(picking, gprim)
//
//    sceneUpdated()
//    updateKeyWindow()

// functions expected by wv
//    wvInitUI()                called by wv-render.js
//    wvUpdateUI()              called by wv-render.js
//    wvServerMessage(text)     called by wv-socket.js
//    wvServerDown()            called by wv-socket.js
//    wvUpdateCanvas(gl)        called by ESP.html
//    wv.sceneUpdated()

// functions associated with button presses (and associated button presses)
//    cmdFile()
//    cmdFileNew()
//    cmdFileOpen()
//    cmdFileExport()
//    cmdFileEdit(e,indx)
//       editorOk()
//       editorCancel()
//    cmdCaps()
//       cmdCapsUpdateIntent()
//       cmdCapsCommitPhase()
//       cmdCapsRevertPhase()
//       cmdCapsListPhases()
//       cmdCapsListAnalyses()
//       cmsCapsListBounds()
//       cmdCapsAttachFile()
//       cmdCapsViewFiles()
//    cmdTool()
//       cmdDone()
//    cmdStepThru(direction)
//    cmdHelp()
//    cmdTest()
//    cmdCollab()
//       cmdCollabListUsers()
//       cmdCollabPass(e,indx)
//       cmdCollabSync()
//       cmdCollabUnsync()
//       cmdCollabMessage()
//    cmdOverlayBeg(timName,overlayName)
//       cmdOverlayEnd()

// functions associated with menu selections (and associated button presses)
//    addCval()
//    editCval(e)
//    addPmtr()
//    editPmtr(e)
//       addRow()
//       addColumn()
//       compGeomSens()
//       compTessSens()
//       setVel(e)
//       clrVels()
//    editValuOk()
//    editValuCancel()
//    delValu()
//    showOutpmtrs()
//       showOutpmtrsOk()
//    addBrch()
//       addBrchOk()
//       addBrchCancel()
//    editBrch(e)
//       addAttr()
//       delBrch(e)
//       delAttr()
//       showBrchAttrs()
//       showBrchArgs()
//       buildTo()
//       editBrchOk()
//       editBrchCancel()
//    allPmtrs()
//    allLocals()
//    allBrchs()
//    chgDisplay()
//    showBodyAttrs(e)

// interface functions for "main" mode
//    main.cmdHome()
//    main.cmdLeft()
//    main.cmdRite()
//    main.cmdBotm()
//    main.cmdTop()
//    main.cmdIn()
//    main.cmdOut()
//
//    main.mouseDown(e)
//    main.mouseMove(e)
//    main.mouseUp(e)
//    main.mouseWheel(e)
//    main.mouseLeftCanvas(e)
//
//    main.keyPress(e)
//    main.keyDown(e)
//    main.keyUp(e)
//
//    main.cmdSolve()
//    main.cmdUndo()

// functions associated with the mouse in the key window
//    setKeyLimits(e)

// functions associated with the mouse in the messages window
//    gotoCsmError(e)

// functions associated with toggling display settings
//    toggleViz(e)
//    toggleGrd(e)
//    toggleTrn(e)
//    toggleOri(e)
//    cancelStepThru()
//    modifyDisplayType(e)
//    modifyDisplayFilter(e)
//    displayFilterOff()

// functions associated with wireframes
//    createWireframes()
//    drawWireframes(gl)

// function associated with lastPoint
//    drawLastPoint(gl)

// functions associated with a Tree in the treefrm
//    Tree(doc, treeId) - constructor
//    TreeAddNode(iparent, name, tooltip, gprim, click,
//                prop1, cbck1,
//                prop2, cbck2,
//                prop3, cbck3)
//    TreeBuild()
//    TreeClear()
//    TreeContract(inode)
//    TreeExpand(inode)
//    TreeProp(inode, iprop, onoff)
//    TreeUpdate()

// helper functions
//    activateBuildButton()
//    resizeFrames()            called by ESP.html
//    browserClosing()          called by ESP.html
//    changeMode(newMode)
//    rebuildTreeWindow(x)
//    checkIfWithBall()
//    checkIfFree()
//    postMessage(mesg)
//    setupEditBrchForm()
//    setupEditValuForm()
//    browserToServer(text)
//    numberOfPmtrChanges()
//    numberOfBrchChanges()
//    unhighlightColumn1()
//    changeProp(inode,icol,state)
//    cmdEditCopy(cm)
//    cmdEditCut(cm)
//    cmdEditPaste(cm)
//    cmdEditFind(cm)
//    cmdEditFindNext(cm)
//    cmdEditFindPrev(cm)
//    cmdEditReplace(cm)
//    cmdEditComment()
//    cmdEditIndent()
//    cmdEditHint(line)
//    cmdEditUndo(cm)
//    cmdEditDebug()
//    printObject(obj)
//    sprintf()
//    CodeMirror.defineSimpleMode(mode, options)

"use strict";


//
// callback when "undoButton" is pressed (called by ESP.html)
//
var cmdUndo = function() {
    if        (wv.overlay !== undefined && wv.overlay.cmdUndo !== undefined) {
        wv.overlay.cmdUndo();
    } else if (wv.curTool.cmdUndo !== undefined) {
        wv.curTool.cmdUndo();
    }
};


//
// callback when "solveButton" is pressed (called by ESP.html)
//
var cmdSolve = function() {

    if (checkIfWithBall() === false) {
        return;
    }

    // if my*Menu is currently posted, delete it/them now
    document.getElementById("myFileMenu"  ).classList.remove("showFileMenu"  );
    document.getElementById("myCapsMenu"  ).classList.remove("showCapsMenu"  );
    document.getElementById("myToolMenu"  ).classList.remove("showToolMenu"  );
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    if        (wv.overlay !== undefined && wv.overlay.cmdSolve !== undefined) {
        wv.overlay.cmdSolve();
    } else if (wv.curTool.cmdSolve !== undefined) {
        wv.curTool.cmdSolve();
    }
};


//
// callback when "saveButton" is pressed (called by ESP.html)
//
var cmdSave = function() {
    if        (wv.overlay !== undefined && wv.overlay.cmdSave !== undefined) {
        wv.overlay.cmdSolve();
    } else if (wv.curTool.cmdSave !== undefined) {
        wv.curTool.cmdSave();
    }
};


//
// callback when "quitButton" is pressed (called by ESP.html)
//
var cmdQuit = function() {
    if        (wv.overlay !== undefined && wv.overlay.cmdQuit !== undefined) {
        wv.overlay.cmdQuit();
    } else if (wv.curTool.cmdQuit !== undefined) {
        wv.curTool.cmdQuit();
    }
};


//
// callback when "homeButton" is pressed (called by ESP.html)
//
var cmdHome = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdHome !== undefined) {
            wv.overlay.cmdHome();
        } else if (wv.curTool.cmdHome !== undefined) {
            wv.curTool.cmdHome();
        } else if (wv.usingMain == 1) {
            main.cmdHome();
        }
    }
};


//
// callback when "leftButton" is pressed (called by ESP.html)
//
var cmdLeft = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdLeft !== undefined) {
            wv.overlay.cmdLeft();
        } else if (wv.curTool.cmdLeft !== undefined) {
            wv.curTool.cmdLeft();
        } else if (wv.usingMain == 1) {
            main.cmdLeft();
        }
    }
};


//
// callback when "riteButton" is pressed (called by ESP.html)
//
var cmdRite = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdRite !== undefined) {
            wv.overlay.cmdRite();
        } else if (wv.curTool.cmdRite !== undefined) {
            wv.curTool.cmdRite();
        } else if (wv.usingMain == 1) {
            main.cmdRite();
        }
    }
};


//
// callback when "botmButton" is pressed (called by ESP.html)
//
var cmdBotm = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdBotm !== undefined) {
            wv.overlay.cmdBotm();
        } else if (wv.curTool.cmdBotm !== undefined) {
            wv.curTool.cmdBotm();
        } else if (wv.usingMain == 1) {
            main.cmdBotm();
        }
    }
};


//
// callback when "topButton" is pressed (called by ESP.html)
//
var cmdTop = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdTop !== undefined) {
            wv.overlay.cmdTop();
        } else if (wv.curTool.cmdTop !== undefined) {
            wv.curTool.cmdTop();
        } else if (wv.usingMain == 1) {
            main.cmdTop();
        }
    }
};


//
// callback when "inButton" is pressed (called by ESP.html)
//
var cmdIn = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdIn !== undefined) {
            wv.overlay.cmdIn();
        } else if (wv.curTool.cmdIn !== undefined) {
            wv.curTool.cmdIn();
        } else if (wv.usingMain == 1) {
            main.cmdIn();
        }
    }
};


//
// callback when "outButton" is pressed (called by ESP.html)
//
var cmdOut = function() {
    if (checkIfFree()) {
        if        (wv.overlay !== undefined && wv.overlay.cmdOut !== undefined) {
            wv.overlay.cmdOut();
        } else if (wv.curTool.cmdOut !== undefined) {
            wv.curTool.cmdOut();
        } else if (wv.usingMain == 1) {
            main.cmdOut();
        }
    }
};


//
// callback when any mouse is pressed in canvas (when wv.usingMain==1)
//
var mouseDown = function(e) {
    if        (wv.overlay !== undefined && wv.overlay.mouseDown !== undefined) {
        wv.overlay.mouseDown(e);
    } else if (wv.curTool.mouseDown !== undefined) {
        wv.curTool.mouseDown(e);
    } else if (wv.usingMain == 1 && wv.myRole != 1) {
        main.mouseDown(e);
    }
};


//
// callback when the mouse moves in canvas (when wv.usingMain==1)
//
var mouseMove = function (e) {
    if        (wv.overlay !== undefined && wv.overlay.mouseMove !== undefined) {
        wv.overlay.mouseMove(e);
    } else if (wv.curTool.mouseMove !== undefined) {
        wv.curTool.mouseMove(e);
    } else if (wv.usingMain == 1 && wv.myRole != 1) {
        main.mouseMove(e);
    }
};


//
// callback when the mouse is released in canvas (when wv.usingMain==1)
//
var mouseUp = function (e) {
    if        (wv.overlay !== undefined && wv.overlay.mouseUp !== undefined) {
        wv.overlay.mouseUp(e);
    } else if (wv.curTool.mouseUp !== undefined) {
        wv.curTool.mouseUp(e);
    } else if (wv.usingMain == 1 && wv.myRole != 1) {
        main.mouseUp(e);
    }
};


//
// callback when the mouse wheel is rolled in canvas (when wv.usingMain==1)
//
var mouseWheel = function (e) {
    if        (wv.overlay !== undefined && wv.overlay.mouseWheel !== undefined) {
        wv.overlay.mouseWheel(e);
    } else if (wv.curTool.mouseWheel !== undefined) {
        wv.curTool.mouseWheel(e);
    } else if (wv.usingMain == 1 && wv.myRole != 1) {
        main.mouseWheel(e);
    }
};


//
// callback when the mouse leaves the canvas (when wv.usingMain==1)
//
var mouseLeftCanvas = function (e) {
    if        (wv.overlay !== undefined && wv.overlay.mouseLeftCanvas !== undefined) {
        wv.overlay.mouseLeftCanvas(e);
    } else if (wv.curTool.mouseLeftCanvas !== undefined) {
        wv.curTool.mouseLeftCanvas(e);
    } else if (wv.usingMain == 1 && wv.myRole != 1) {
        main.mouseLeftCanvas(e);
    }
};


//
// callback when a key is pressed (when wv.usingMain==1)
//
var keyPress = function (e) {
    if (wv.curTool == main && wv.overlay === undefined) {
        main.keyPress(e);
    } else if (wv.overlay !== undefined && wv.overlay.keyPress !== undefined) {
        var handled = wv.overlay.keyPress(e);

        if (handled == 0) {
            main.keyPress(e);
        }
    } else if (wv.curTool.keyPress !== undefined) {
        var handled = wv.curTool.keyPress(e);

        if (handled == 0) {
            main.keyPress(e);
        }
    } else if (wv.usingMain == 1) {
        main.keyPress(e);
    }
};


//
// callback when an arrow... or shift is pressed (needed for Chrome)
//
var keyDown = function (e) {
    if        (wv.overlay !== undefined && wv.overlay.keyDown !== undefined) {
        wv.overlay.keyDown(e);
    } else if (wv.curTool.keyDown != undefined) {
        wv.curTool.keyDown(e);
    } else if (wv.usingMain == 1) {
        main.keyDown(e);
    }
};


//
// callback when a shift is released (needed for Chrome)
//
var keyUp = function (e) {
    if        (wv.overlay !== undefined && wv.overlay.keyUp !== undefined) {
        wv.overlay.keyUp(e);
    } else if (wv.curTool.keyUp !== undefined) {
        wv.curTool.keyUp(e);
    } else if (wv.usingMain == 1) {
        main.keyUp(e);
    }
};

////////////////////////////////////////////////////////////////////////////////


//
// callback when the user interface is to be initialized (called by wv-render.js)
//
var wvInitUI = function () {
    // alert("in wvInitUI()");

    // set up extra storage for matrix-matrix multiplies
    wv.uiMatrix   = new J3DIMatrix4();
    wv.saveMatrix = new J3DIMatrix4(wv.mvMatrix);   // matrix used in save operations
    wv.wireMatrix = new J3DIMatrix4();              // matrix used in wireframe operations

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
    wv.buildTree =  0;             // =1 if tree is being rebuilt
    wv.bodynames = undefined;      // list of current Body names
    wv.wireframe = undefined;      // corners of each wireframe
    wv.debugUI   =  0;             // set to 1 for console messages
    wv.idntStat  =  0;             // -1 server is identified
                                   //  0 need to identify server
                                   // >0 waiting for server to identify
    wv.cvalStat  =  0;             // -2 latest Caps Values are in Tree
                                   // -1 latest Caps Values not in Tree (yet)
                                   //  0 need to request CapsValues
                                   // >0 waiting for Caps Values (request already made)
    wv.pmtrStat  =  0;             // -2 latest Parameters are in Tree
                                   // -1 latest Parameters not in Tree (yet)
                                   //  0 need to request Parameters
                                   // >0 waiting for Parameters (request already made)
    wv.brchStat  =  0;             // -2 latest Branches are in Tree
                                   // -1 latest Branches not in Tree (yet)
                                   //  0 need to request Branches
                                   // >0 waiting for Branches (request already made)
    wv.builtTo   = 99999;          // last Branch in previous successful build
    wv.menuEvent = undefined;      // event associated with click in Tree
    wv.server    = undefined;      // string passed back from "identify;"
    wv.myRole    = 0;              // =0 hasBall, =1 sync'd, =2 unsync'd
    wv.numUsers  = 1;              // number of users
    wv.userNames = "";             // bar-separatede list of userNames
    wv.myName    = "*host*";       // my username
    wv.lastXform = null;           // last xfrom received while not sync'd
    wv.plotType  =  0;             // =0 mono, =1 ubar, =2 vbar, =3 cmin, =4 cmax, =5 gc, =6 normals, =10 erep, =11 plugs, =12 pyscript, =13 mitten
    wv.loLimit   = -1;             // lower limit in key
    wv.upLimit   = +1;             // upper limit in key
    wv.nchanges  = 0;              // number of Branch or Parameter changes by browser
    wv.pyname    = "";             // name of the .py file that started ESP
    wv.filenames = "|";            // name of the .csm file (and .udc files)
    wv.fileintro = false;          // =true if intro has been posted
    wv.fileindx  = undefined;      // index of file being editted: -1 <new file>, 0 *.csm, >0 *.udc
    wv.linenum   = 0;              // line number to start editing
    wv.lastfile  = "";             // last file that was editted
    wv.lastline  = -1;             // line number in last file editted
    wv.curFile   = "";             // contents of current .csm file
    wv.codeMirror= undefined;      // pointer to CodeMirror object
    wv.clipboard = "";             // clipboard for copy/cut/paste
    wv.curMode   =  0;             //-1 to disable buttons, etc.
                                   // 0 show WebViewer in canvas
                                   // 1 show addBrchForm
                                   // 2 show editBrchForm with addBrchHeader
                                   // 3 show editBrchForm with editBrchHeader
                                   // 4 show editValuForm with addValuHeader
                                   // 5 show editValuForm with editValuHeader
                                   // 6 show showOutpmtrsForm
                                   // 7 show editorForm
                                   // 8 show sketcherForm and run sketch
                                   // 9 show glovesForm
                                   // 10 show WebViewer in canvas and run ereped
                                   // 11 show WebViewer in canvas and run plugs
                                   // 12 show WebViewer in canvas and run pyscript
                                   // 13 show WebViewer in canvas and run mitten
                                   // 14 show sketcherForm and run plotter
    wv.curTool   = main;           // current tool
    wv.timName   = "";             // name of TIM being held
    wv.overlay   = undefined;      // name of process that overlays the current mode
    wv.usingMain =  1;             // =1 if using 3D graphics window
    wv.capsIntent = undefined;     // intent of current CAPS session
    wv.curStep   =  0;             // >0 if in StepThru mode
    wv.curCval   = -1;             // Caps Value being editted (or -1)
    wv.curPmtr   = -1;             // Parameter being editted (or -1)
    wv.curBrch   = -1;             // Branch being editted (or -1)
    wv.afterBrch = -1;             // Branch to add after (or -1)
    wv.numArgs   =  0;             // number of arguments in wv.curBrch
    wv.sgData    = {};             // scene graph metadata
    wv.scale     =  1;             // scale factor for axes
    wv.getFocus  = undefined;      // entry to get first focus
    wv.lastPoint = undefined;      // array of xyz at last @ command
//  wv.centerV                     // set to 1 to center view and rotation
//  wv.pick                        // set to 1 to turn picking on
//  wv.locate                      // set to 1 to turn locating on
//  wv.sceneGraph                  // pointer to sceneGraph
//  wv.sceneUpdated()              // callback when scene graphs has been updated
//  wv.picked                      // sceneGraph object that was picked
//  wv.located                     // sceneGraph object that was located
//  wv.sceneUpd                    // should be set to 1 to re-render scene
//  wv.sgUpdate                    // =1 if the sceneGraph has been updated
//  wv.canvasID                    // ID of main canvas
//  wv.canvasKY                    // ID of key  canvas
//  wv.drawKey                     // =1 if key window needs to be redrawn
//  wv.socketUt.send(text)         // function to send text to server
//  wv.plotAttrs                   // plot attributes

    window.addEventListener('beforeunload',  browserClosing,  false);

    document.addEventListener('keypress',   keyPress,         false);
    document.addEventListener('keydown',    keyDown,          false);
    document.addEventListener('keyup',      keyUp,            false);

    var canvas = document.getElementById(wv.canvasID);
    canvas.addEventListener(  'mousedown',  mouseDown,        false);
    canvas.addEventListener(  'mousemove',  mouseMove,        false);
    canvas.addEventListener(  'mouseup',    mouseUp,          false);
    canvas.addEventListener(  "wheel",      mouseWheel,       false);
    canvas.addEventListener(  'mouseout',   mouseLeftCanvas,  false);

    var sketcher = document.getElementById("sketcher");
    sketcher.addEventListener('mousedown',  sketch.mouseDown, false);
    sketcher.addEventListener('mousemove',  sketch.mouseMove, false);
    sketcher.addEventListener('mouseup',    sketch.mouseUp,   false);

    var gloves = document.getElementById("gloves");
    gloves.addEventListener(  'mousedown',  gloves.mouseDown, false);
    gloves.addEventListener(  'mousemove',  gloves.mouseMove, false);
    gloves.addEventListener(  'mouseup',    gloves.mouseUp,   false);

    var keycan = document.getElementById(wv.canvasKY);
    keycan.addEventListener(  'mouseup',    setKeyLimits,     false);

    var msgwin = document.getElementById("brframe");
    msgwin.addEventListener(  'dblclick',   gotoCsmError,     false);

    document.getElementById("exitOlayBtn").hidden = true;
};


//
// callback when the user interface should be updated (called by wv-render.js)
//
var wvUpdateUI = function () {
    // alert("in wvUpdateUI()");

    // special code if ESP was started with a .py file
    if (wv.pyname != "") {
        var temp = wv.pyname;
        wv.pyname = "";

        pyscript.launch(temp);

        return;
    }

    // special code for delayed-picking mode
    if (wv.picking > 0) {

        // if something is picked, post a message
        if (wv.picked !== undefined) {

            // second part of 'A' operation
            if (wv.picking == 65) {
                while (1) {             // used to jump out on error

                    // get Attribute name and value
                    var newAttrName = prompt("Enter new Attribute name:");
                    if (newAttrName === null) {
                        break;
                    }

                    var newAttrValue = prompt("Enter new Attribute value" +
                                              "(either a semi-colon separated list" +
                                              " or begin with a '$' for a string)");
                    if (newAttrValue === null) {
                        break;
                    }

                    // add statements after last executed Branch or at end of Branches
                    if (wv.builtTo < 99999) {
                        var ibrch = Number(wv.builtTo);
                    } else {
                        var ibrch = brch.length;
                    }

                    // determine the Body to select
                    var ibody = -1;
                    var gprimList = wv.picked.gprim.trim().split(" ");
                    if (gprimList.length == 3) {
                        var tempList = wv.sgData[gprimList[0]];
                        for (var ii = 0; ii < tempList.length; ii+=2) {
                            if (tempList[ii] == "_body") {
                                ibody = Number(tempList[ii+1]);
                                break;
                            }
                        }
                    } else if (gprimList.length == 4) {
                        ibody = Number(gprimList[1]);
                    }
                    if (ibody <= 0) {
                        alert("Could not find the Body");
                        break;
                    }

                    // find the _faceID, _edgeID, or _nodeID Attributes
                    var attrs = wv.sgData[wv.picked.gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrs[i] == "_faceID") {

                            // create the "select body" statement
                            browserToServer("newBrch|"+ibrch+"|select|$body|"+ibody+"|||||||||");

                            // create the "select face" statement
                            var faceIDlist = attrs[i+1].trim().split(" ");
                            browserToServer("newBrch|"+(ibrch+1)+"|select|$face|"+faceIDlist[0]+"|"+
                                                                                  faceIDlist[1]+"|"+
                                                                                  faceIDlist[2]+"||||||||");

                            // add the Attribute to the new select statement
                            browserToServer("setAttr|"+(ibrch+2)+"|"+newAttrName+"|1|"+newAttrValue+"|");

                            postMessage("Adding: select $body "+ibody+"\n" +
                                        "        select $face "+faceIDlist[0]+" "+
                                                                faceIDlist[1]+" "+
                                                                faceIDlist[2]+"\n" +
                                        "            attribute "+newAttrName+" "+newAttrValue+"\n" +
                                        "====> Re-build is needed <====");

                            // get an updated version of the Branches and activate Build button
                            wv.brchStat = 0;

                            activateBuildButton();
                            break;
                        } else if (attrs[i] == "_edgeID") {

                            // create the "select body" statement
                            browserToServer("newBrch|"+ibrch+"|select|$body|"+ibody+"|||||||||");

                            // create the "select edge" statement
                            var edgeIDlist = attrs[i+1].trim().split(" ");
                            browserToServer("newBrch|"+(ibrch+1)+"|select|$edge|"+edgeIDlist[0]+"|"+
                                                                                  edgeIDlist[1]+"|"+
                                                                                  edgeIDlist[2]+"|"+
                                                                                  edgeIDlist[3]+"|"+
                                                                                  edgeIDlist[4]+"||||||");

                            // add the Attribute to the new select statement
                            browserToServer("setAttr|"+(ibrch+2)+"|"+newAttrName+"|1|"+newAttrValue+"|");

                            postMessage("Adding: select $body "+ibody+"\n" +
                                        "        select $edge "+edgeIDlist[0]+" "+
                                                                edgeIDlist[1]+" "+
                                                                edgeIDlist[2]+" "+
                                                                edgeIDlist[3]+" "+
                                                                edgeIDlist[4]+"\n" +
                                        "            attribute "+newAttrName+" "+newAttrValue+"\n" +
                                        "====> Re-build is needed <====");

                            // get an updated version of the Branches and activate Build button
                            wv.brchStat = 0;

                            activateBuildButton();
                            break;
                        } else if (attrs[i] == "_nodeID") {

                            // create the "select body" statement
                            browserToServer("newBrch|"+ibrch+"|select|$body|"+ibody+"|||||||||");

                            // create the "select edge" statement
                            var edgeIDlist = attrs[i+1].trim().split(" ");
                            browserToServer("newBrch|"+(ibrch+1)+"|select|$node|"+attrs[i+1]+"||||||||||");

                            // add the Attribute to the new select statement
                            browserToServer("setAttr|"+(ibrch+2)+"|"+newAttrName+"|1|"+newAttrValue+"|");

                            postMessage("Adding: select $body "+ibody+"\n" +
                                        "        select $node "+attrs[i+1]+"\n"+
                                        "            attribute "+newAttrName+" "+newAttrValue+"\n" +
                                        "====> Re-build is needed <====");

                            // get an updated version of the Branches and activate Build button
                            wv.brchStat = 0;

                            activateBuildButton();
                            break;
                        }
                    }
                    break;
                }

            // second part of 'g' operation
            } else if (wv.picking == 103) {
                postMessage("Toggling grid of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
                            changeProp(inode, 2, "on");
                        } else {
                            changeProp(inode, 2, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of 'o' operation
            } else if (wv.picking == 111) {
                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if (wv.picked.gprim.search(" Edge ") >= 0) {
                            postMessage("Toggling orientation of "+wv.picked.gprim);

                            if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
                                changeProp(inode, 3, "on");
                            } else {
                                changeProp(inode, 3, "off");
                            }
                            myTree.update();
                            break;
                        } else {
                            postMessage("Cannot toggle orientation of a non-Edge");
                        }
                    }
                }

            // second part of 't' operation
            } else if (wv.picking == 116) {
                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if (wv.picked.gprim.search(" Face ") >= 0) {
                            postMessage("Toggling transparency of "+wv.picked.gprim);

                            if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
                                changeProp(inode, 3, "on");
                            } else {
                                changeProp(inode, 3, "off");
                            }
                            myTree.update();
                            break;
                        } else {
                            postMessage("Cannot toggle transparency of a non-Face");
                        }
                    }
                }

            // second part of 'v' operation
            } else if (wv.picking == 118) {
                postMessage("Toggling visibility of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
                            changeProp(inode, 1, "on");
                        } else {
                            changeProp(inode, 1, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of '^' operation
            } else if (wv.picking == 94) {
                var mesg = "Picked: "+wv.picked.gprim;

                try {
                    var attrs = wv.sgData[wv.picked.gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        mesg = mesg + "\n        "+attrs[i]+"= "+attrs[i+1];
                    }
                } catch (x) {
                }
                postMessage(mesg);

            // second part of another operation (to be processed by wv.curTool)
            } else if (wv.curTool.keyPressPart2 !== undefined) {
                wv.curTool.keyPressPart2(wv.picking, wv.picked.gprim);
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

                if (wv.lastPoint === undefined) {
                    postMessage("Approximate location: x="+xloc.toFixed(4)+", y="+yloc.toFixed(4)+
                                ", z="+zloc.toFixed(4));
                } else {
                    var dist = Math.sqrt(Math.pow(xloc-wv.lastPoint[0], 2)
                                       + Math.pow(yloc-wv.lastPoint[1], 2)
                                       + Math.pow(zloc-wv.lastPoint[2], 2));
                    postMessage("Approximate location: x="+xloc.toFixed(4)+", y="+yloc.toFixed(4)+
                                ", z="+zloc.toFixed(4)+", dist="+dist.toFixed(6));
                }

                if (wv.lastPoint !== undefined            &&
                    Math.abs(xloc-wv.lastPoint[0]) < 1e-6 &&
                    Math.abs(yloc-wv.lastPoint[1]) < 1e-6 &&
                    Math.abs(zloc-wv.lastPoint[2]) < 1e-6   ) {
                    wv.lastPoint = undefined;
                } else {
                    wv.lastPoint = [xloc, yloc, zloc];
                }

                if (wv.myRole == 0 && wv.numUsers > 1) {
                    if (wv.lastPoint !== undefined) {
                        browserToServer("lastPoint|"+wv.lastPoint[0]+"|"+wv.lastPoint[1]+"|"+wv.lastPoint[2]+"|");
                    } else {
                        browserToServer("lastPoint|off|");
                    }
                }
                wv.sceneUpd  = 1;
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

    // if the server is not identified, try to identify it now
    if (wv.idntStat > 0) {
        wv.idntStat--;
    } else if (wv.idntStat == 0) {
        try {
            browserToServer("identify|");
            wv.idntStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.idntStat = 10;
       }
    }

    // if the Caps Values are scheduled to be updated, send a message to
    //    get the Caps Values now
    if (wv.cvalStat > 0) {
        wv.cvalStat--;
    } else if (wv.server === undefined) {
    } else if (wv.server != "serveCAPS") {
        wv.cvalStat = -2;
    } else if (wv.cvalStat == 0) {
        try {
            browserToServer("getCvals|1|");
            wv.cvalStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.cvalStat = 10;
        }
    }

    // if the Parameters are scheduled to be updated, send a message to
    //    get the Parameters now
    if (wv.pmtrStat > 0) {
        wv.pmtrStat--;
    } else if (wv.pmtrStat == 0) {
        try {
            browserToServer("getPmtrs|");
            wv.pmtrStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.pmtrStat = 10;
        }
    }

    // if the Branches are scheduled to be updated, send a message to
    //    get the Branches now
    if (wv.brchStat > 0) {
        wv.brchStat--;
    } else if (wv.server === undefined) {
    } else if (wv.server == "serveCAPS") {
        wv.brchStat = -2;
    } else if (wv.brchStat == 0) {
        try {
            browserToServer("getBrchs|1|");
            wv.brchStat = -1;
        } catch (e) {
            // could not send command, so try again after 10 cycles
            wv.brchStat = 10;
        }
    }

    // if the scene graph has been updated, (re-)build the Tree
    if (wv.sgUpdate == 1) {

        if (wv.sceneGraph === undefined) {
            alert("wv.sceneGraph is undefined --- but we need it");
        }

        rebuildTreeWindow(1);

    // if the Caps Values, Parameters, or Branches have been updated, (re-) build the Tree
    } else if (wv.cvalStat == -1 || wv.pmtrStat == -1 || wv.brchStat == -1) {

        if (wv.sceneGraph === undefined) {
            alert("wv.sceneGraph is undefined --- but we need it");
        }

        rebuildTreeWindow(2);
    }

    // deal with key presses
    if (wv.keyPress != -1 && wv.usingMain == 1) {

        var myKeyPress = String.fromCharCode(wv.keyPress);

        // check to see if wv.curTool can process this
        var done = 0;

        if (wv.overlay !== undefined && wv.overlay.keyPress !== undefined) {
            done = wv.overlay.keyPress(myKeyPress);
        } else if (wv.curTool.keyPressPart1 !== undefined) {
            done = wv.curTool.keyPressPart1(myKeyPress);
        }

        if (done == 1) {

        // '?' -- help
        } else if (myKeyPress == "?") {
            postMessage("........................... Viewer Cursor options ...........................\n" +
                        "ctrl-h <Home> - initial view             ctrl-f        - front view          \n" +
                        "ctrl-l        - leftside view            ctrl-r        - riteside view       \n" +
                        "ctrl-t        - top view                 ctrl-b        - bottom view         \n" +
                        "ctrl-i <PgUp> - zoom in                  ctrl-o <PgDn> - zoom out            \n" +
                        "<Left>        - rotate or xlate left     <Rite>        - rotate or xlate rite\n" +
                        "<Up>          - rotate or xlate up       <Down>        - rotate or xlate down\n" +
                        ">             - save view                <             - recall view         \n" +
                        "ctrl-> .      - save view to file        ctrl-< ,      - read view from file \n" +
                        "^ 6           - query object at cursor   @ 2           - get coords. @ cursor\n" +
                        "v             - toggle Viz at cursor     g             - toggle Grd at cursor\n" +
                        "t             - toggle Trn at cursor     o             - toggle Ori at cursor\n" +
                        "A             - add Attribute at cursor  * 8           - center view @ cursor\n" +
                        "!             - toggle flying mode       ?             - get help            \n" +
                        ".............................................................................");

        // 'A' -- add Attribute at cursor
        } else if (myKeyPress == 'A' && checkIfWithBall()) {
            wv.picking  = 65;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'g' -- toggle grid at cursor
        } else if (myKeyPress == 'g' && wv.modifier == 0 && checkIfFree()) {
            wv.picking  = 103;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'o' -- orientation at cursor
        } else if (myKeyPress == 'o' && wv.modifier == 0 && checkIfFree()) {
            wv.picking  = 111;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 't' -- toggle transparency at cursor
        } else if (myKeyPress == 't' && wv.modifier == 0 && checkIfFree()) {
            wv.picking  = 116;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'v' -- toggle visibility at cursor
        } else if (myKeyPress == 'v' && wv.modifier == 0 && checkIfFree()) {
            wv.picking  = 118;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'n' -- NextStep if in StepThru mode
        } else if (myKeyPress == 'n' && wv.modifier == 0 && checkIfWithBall()) {
            if (wv.curStep > 0) {
                cmdStepThru(+1);
            } else {
                postMessage("Press \"StepThru\" to enter StepThru mode");
            }

        // 'p' -- PrevStep if in StepThru mode
        } else if (myKeyPress == 'p' && wv.modifier == 0 && checkIfWithBall()) {
            if (wv.curStep > 0) {
                cmdStepThru(-1);
            } else {
                postMessage("Press \"StepThru\" to enter StepThru mode");
            }

        // 'f' -- FirstStep if in StepThru mode
        } else if (myKeyPress == 'f' && wv.modifier == 0 && checkIfWithBall()) {
            if (wv.curStep > 0) {
                cmdStepThru(-2);
            } else {
                postMessage("Press \"StepThru\" to enter StepThru mode");
            }

        // 'l' -- LastStep if in StepThru mode
        } else if (myKeyPress == 'l' && wv.modifier == 0 && checkIfWithBall()) {
            if (wv.curStep > 0) {
                cmdStepThru(+2);
            } else {
                postMessage("Press \"StepThru\" to enter StepThru mode");
            }

        // 'x' - bump point locations in wireframe
        } else if ((myKeyPress == 'x' || myKeyPress == 'X' ||
                    myKeyPress == 'y' || myKeyPress == 'Y' ||
                    myKeyPress == 'z' || myKeyPress == 'Z'   ) && checkIfWithBall()) {

            if (wv.wireframe === undefined || wv.wireMatrix === undefined) {
                postMessage("\""+myKeyPress+"\" is only valid with wireframes");
            } else {
//                postMessage("wv.cursorX    "+wv.cursorX);
//                postMessage("wv.cursorY    "+wv.cursorY);
//                postMessage("wv.width      "+wv.width);
//                postMessage("wv.height     "+wv.height);
//                postMessage("wv.focus      "+wv.focus[0]+" "+wv.focus[1]+" "+wv.focus[2]+" "+wv.focus[3])

                // get location of key press (between -1 and +1)
                var scrX = 2.0 * wv.cursorX / (wv.width  - 1.0) - 1.0;
                var scrY = 2.0 * wv.cursorY / (wv.height - 1.0) - 1.0;

                // find closest point
                var dxybest = 9999999;
                var ibest   = -1;
                var myMatrix = wv.wireMatrix.getAsArray();

//                postMessage("mvpMatrix  "+myMatrix[ 0]+" "+myMatrix[ 4]+" "+myMatrix[ 8]+" "+myMatrix[12]);
//                postMessage("           "+myMatrix[ 1]+" "+myMatrix[ 5]+" "+myMatrix[ 9]+" "+myMatrix[13]);
//                postMessage("           "+myMatrix[ 2]+" "+myMatrix[ 6]+" "+myMatrix[10]+" "+myMatrix[14]);
//                postMessage("           "+myMatrix[ 3]+" "+myMatrix[ 7]+" "+myMatrix[11]+" "+myMatrix[15]);

                for (var i = 0; i < wv.wireframe.length/3; i++) {
                    var Xtemp = (wv.wireframe[3*i  ] - wv.focus[0]) / wv.focus[3];
                    var Ytemp = (wv.wireframe[3*i+1] - wv.focus[1]) / wv.focus[3];
                    var Ztemp = (wv.wireframe[3*i+2] - wv.focus[2]) / wv.focus[3];

                    var Xscr = ( myMatrix[ 0] * Xtemp + myMatrix[ 4] * Ytemp + myMatrix[ 8] * Ztemp + myMatrix[12]);
                    var Yscr = ( myMatrix[ 1] * Xtemp + myMatrix[ 5] * Ytemp + myMatrix[ 9] * Ztemp + myMatrix[13]);
                    var Wscr = ( myMatrix[ 3] * Xtemp + myMatrix[ 7] * Ytemp + myMatrix[11] * Ztemp + myMatrix[15]);

                    var dx = scrX - Xscr / Wscr;
                    var dy = scrY - Yscr / Wscr;

//                    postMessage(wv.wireframe[3*i]+" "+wv.wireframe[3*i+1]+" "+wv.wireframe[3*i+2]+"->"+(Xscr/Wscr+1)*(wv.width-1)/2+" "+(Yscr/Wscr+1)*(wv.height-1)/2);

                    var dxytest = dx * dx + dy * dy;
                    if (dxytest < dxybest) {
                        dxybest = dxytest;
                        ibest   = i;
                    }
                }

                if (ibest >= 0) {
                    if        (myKeyPress == 'x') {
                        wv.wireframe[3*ibest  ] -= 0.1;
                        wv.sceneUpd = 1;
                    } else if (myKeyPress == 'X') {
                        wv.wireframe[3*ibest  ] += 0.1;
                        wv.sceneUpd = 1;
                    } else if (myKeyPress == 'y') {
                        wv.wireframe[3*ibest+1] -= 0.1;
                        wv.sceneUpd = 1;
                    } else if (myKeyPress == 'Y') {
                        wv.wireframe[3*ibest+1] += 0.1;
                        wv.sceneUpd = 1;
                    } else if (myKeyPress == 'z') {
                        wv.wireframe[3*ibest+2] -= 0.1;
                        wv.sceneUpd = 1;
                    } else if (myKeyPress == 'Z') {
                        wv.wireframe[3*ibest+2] += 0.1;
                        wv.sceneUpd = 1;
                    }
                }
            }

        // '^' or '6' -- query at cursor
        } else if ((myKeyPress == '^' || myKeyPress == '6') && checkIfFree()) {
            wv.picking  = 94;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '@' or '2' -- locate at cursor
        } else if ((myKeyPress == '@' || myKeyPress == '2') && checkIfFree()) {
            wv.locating = 64;
            wv.locate   = 1;

        // '*' or '8' -- center view
        } else if ((myKeyPress == '*' || myKeyPress == '8') && checkIfFree()) {
            wv.centerV = 1;
            postMessage("View being centered");

        // '!' -- toggle flying mode
        } else if (myKeyPress == "!" && checkIfFree()) {
            if (wv.flying <= 1) {
                postMessage("Turning flying mode ON");
                wv.flying = 10;
            } else {
                postMessage("Turning flying mode OFF");
                wv.flying = 1;
            }

        // '>' -- save view
        } else if (myKeyPress == ">" && checkIfFree()) {
            postMessage("Saving current view");
            wv.saveMatrix.load(wv.mvMatrix);
            wv.sceneUpd = 1;

        // '<' -- recall view
        } else if (myKeyPress == "<" && checkIfFree()) {
            postMessage("Restoring saved view");
            wv.mvMatrix.load(wv.saveMatrix);
            wv.sceneUpd = 1;

        // C-'>' or '.' -- save view to file
        } else if (((wv.keyPress == 46  && wv.modifier == 5) ||
                    (myKeyPress  == "." && wv.modifier == 0)   ) && checkIfWithBall()) {
            var filename = prompt("Enter view filename to save:", "save.view");
            if (filename !== null) {
                postMessage("Saving view to \"" + filename + "\"");
                browserToServer("saveView|"+filename+"|"+wv.scale+"|"+wv.mvMatrix.getAsArray()+"|");
            }

        // C-'<' or ',' -- read view from file
        } else if (((wv.keyPress == 44  && wv.modifier == 5) ||
                    (myKeyPress  == "," && wv.modifier == 0)   ) && checkIfWithBall()) {
            var filename = prompt("Enter view filename to read:", "save.view");
            if (filename !== null) {
                postMessage("Reading view from \"" + filename + "\"");
                browserToServer("readView|"+filename+"|");
            }

        // '<esc>' -- not used
//      } else if (wv.keyPress == 0 && wv.keyCode == 27) {
//          postMessage("<Esc> is not supported.  Use '?' for help");

        // '<Home>' -- initial view
        } else if (wv.keyPress == 0 && wv.keyCode == 36 && checkIfFree()) {
            wv.mvMatrix.makeIdentity();
            wv.scale    = 1;
            wv.sceneUpd = 1;

        // '<end>' -- not used
        } else if (wv.keyPress == 0 && wv.keyCode == 35) {
            postMessage("<End> is not supported.  Use '?' for help");

        // '<PgUp>' -- zoom in
        } else if (wv.keyPress == 0 && wv.keyCode == 33 && checkIfFree()) {
            if (wv.modifier == 0) {
                wv.mvMatrix.scale(2.0, 2.0, 2.0);
                wv.scale *= 2.0;
            } else {
                wv.mvMatrix.scale(1.25, 1.25, 1.25);
                wv.scale *= 1.25;
            }
            wv.sceneUpd = 1;

        // '<PgDn>' -- zoom out
        } else if (wv.keyPress == 0 && wv.keyCode == 34 && checkIfFree()) {
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
        } else if (wv.keyPress == 0 && wv.keyCode == 37 && checkIfFree()) {
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
        } else if (wv.keyPress == 0 && wv.keyCode == 39 && checkIfFree()) {
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
        } else if (wv.keyPress == 0 && wv.keyCode == 38 && checkIfFree()) {
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
        } else if (wv.keyPress == 0 && wv.keyCode == 40 && checkIfFree()) {
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
        } else if (((wv.keyPress == 104 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==   8 && wv.modifier == 4 && wv.keyCode ==  8)   ) && checkIfFree()) {
            main.cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-i' - zoom in (same as <PgUp> without shift)
        } else if (((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) && checkIfFree()) {
            main.cmdIn();
            wv.keyPress = -1;
            return;

        // '+' - zoom in (same as <PgUp> without shift)
        } else if (wv.keyPress ==  43 && wv.modifier == 1 && checkIfFree()) {
            main.cmdIn();
            wv.keyPress = -1;
            return;

        // 'ctrl-o' - zoom out (same as <PgDn> without shift)
        } else if (((wv.keyPress == 111 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==  15 && wv.modifier == 4 && wv.keyCode == 15)   ) && checkIfFree()) {
            main.cmdOut();
            wv.keyPress = -1;
            return;

        // '-' - zoom out (same as <PgDn> without shift)
        } else if (wv.keyPress ==  45 && wv.modifier == 0 && checkIfFree()) {
            main.cmdOut();
            wv.keyPress = -1;
            return;

        // 'ctrl-f' - front view (same as <Home>)
        } else if (((wv.keyPress == 102 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==   6 && wv.modifier == 4 && wv.keyCode ==  6)   ) && checkIfFree()) {
            main.cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-r' - riteside view
        } else if (((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) && checkIfFree()) {
            main.cmdRite();
            wv.keyPress = -1;
            return;

        // 'ctrl-l' - leftside view
        } else if (((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) && checkIfFree()) {
            main.cmdLeft();
            wv.keyPress = -1;
            return;

        // 'ctrl-t' - top view
        } else if (((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) && checkIfFree()) {
            main.cmdTop();
            wv.keyPress = -1;
            return;

        // 'ctrl-b' - bottom view
        } else if (((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
                    (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) && checkIfFree()) {
            main.cmdBotm();
            wv.keyPress = -1;
            return;

        // NOP
        } else if (wv.keyPress == 0 && wv.modifier == 0) {

//$$$        // unknown command
//$$$      } else {
//$$$          postMessage("Unknown command (keyPress="+wv.keyPress
//$$$                      +", modifier="+wv.modifier
//$$$                      +", keyCode="+wv.keyCode+").  Use '?' for help");
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
};


//
// callback when a (text) message is received from the server (called by wv-socket.js)
//
var wvServerMessage = function (text) {
    // alert("in wvServerMessage(text="+text+")");

    // remove trailing NULL
    if (text.charCodeAt(text.length-1) == 0) {
        text = text.substring(0, text.length-1);
    }

    if (wv.debugUI) {
        var date = new Date;
        console.log("("+date.toTimeString().substring(0,8)+") browser<--server: "+text.substring(0,40));
    }

    // if it is a null message, do nothing
    if (text.length < 1) {

    // if it starts with "identify|" post a message */
    } else if (text.substring(0,9) == "identify|") {
        if (wv.server === undefined) {
            var textList = text.split("|");

            wv.server = textList[1];
            postMessage("ESP has been initialized and is attached to '"+wv.server+"'");

            /* get the Cvals, Pmtrs, and Brchs */
            if (wv.server == "serveCAPS") {
                browserToServer("getCvals|0|");
                wv.cvalStat = -1;
            }

            browserToServer("getPmtrs|");
            wv.pmtrStat = -1;

            if (wv.server != "serveCAPS") {
                browserToServer("getBrchs|");
                wv.brchStat = -1;
            }

            if (wv.server != "serveCSM" && wv.server != "serveESP" && wv.server != "serveCAPS") {
                alert("You must be attached to \"serveCSM\" or \"serveESP\" or \"serveCAPS\"");
            }

            if (Number(textList[2]) > 1) {
                wv.myName = null;
                while (wv.myName === null) {
                    wv.myName = prompt("You are connection "+textList[2]+" to the server\n" +
                                       "Enter your username:");
                }
            }
            browserToServer("userName|"+wv.myName+"|");

            if (textList[3].length > 0) {
                wv.pyname = textList[3];
            }

            if (wv.server == "serveCAPS") {
                wv.capsProj  = textList[4];
                wv.capsPhase = textList[5];

                wv.capsIntent = prompt("Working on Phase "+wv.capsPhase+" of Project \""+wv.capsProj+"\"\n"+
                                       "What is the intent of the Phase?");
            }

            if (wv.server == "serveCAPS") {
                document.getElementById("fileMenuBtn").hidden = true;
                document.getElementById("capsMenuBtn").hidden = false;
            } else {
                document.getElementById("fileMenuBtn").hidden = false;
                document.getElementById("capsMenuBtn").hidden = true;
            }
        }

    // if it starts with "userName|" store the list of current users
    } else if (text.substring(0,9) == "userName|") {
        var textList    = text.split("|");
        var oldTextList = wv.userNames.split("|");

        // set up the userNames
        wv.userNames = "|";
        wv.numUsers  = 0;
        for (var i = 2; i < textList.length-1; i++) {
            wv.userNames += textList[i]+"|";
            wv.numUsers++;
        }

        // inform remaining users if someone has left the collaboration
        for (var iold = 1; iold < oldTextList.length-1; iold++) {
            var okay = 0;
            for (var inew = 2; inew < textList.length-1; inew++) {
                if (oldTextList[iold] == textList[inew]) {
                    okay = 1;
                }
            }
            if (okay == 0) {
                postMessage(oldTextList[iold]+" has left the collaboration");
            }
        }

        // inform if someone has joined the collaboration
        if (oldTextList.length > 1) {
            for (var inew = 2; inew < textList.length-1; inew++) {
                okay = 0;
                for (var iold = 1; iold < oldTextList.length-1; iold++) {
                    if (oldTextList[iold] == textList[inew]) {
                        okay = 1;
                    }
                }
                if (okay == 0) {
                    postMessage(textList[inew]+" has joined the collaboration");
                }
            }
        }

        // I am the only user
        if (textList.length == 4) {
            wv.myRole = 0;

            var button = document.getElementById("collabMenuBtn");
            button.hidden = true;

            try {
                wv.codeMirror.setOption("readOnly", false);
            } catch (x) {
            }

        // there are multiple users, but I have the ball
        } else if (wv.myName == textList[Number(textList[1])+2]) {
            wv.myRole = 0;

            var button = document.getElementById("collabMenuBtn");
            button.hidden = false;
            button.style.backgroundColor = "#AFFFAF";    // greenish

            try {
                wv.codeMirror.setOption("readOnly", false);
            } catch (x) {
            }

        // there are multiple users and I do not have the ball
        } else {
            postMessage("\""+textList[Number(textList[1])+2]+"\" has the ball");

            // I previously had the ball, so revert to un-sync'd
            if (wv.myRole == 0) {
                wv.myRole = 2;
            }

            var button = document.getElementById("collabMenuBtn");
            button.hidden = false;

            if (wv.myRole == 1) {
                button.style.backgroundColor = "#FFFFAF";    // yellow
            } else {
                button.style.backgroundColor = null;         // white
            }

            try {
                wv.codeMirror.setOption("readOnly",  true);
            } catch (x) {
            }
        }

        browserToServer("getFilenames|");

    // if it starts with "sgData|" store the auxiliary scene graph data
    } else if (text.substring(0,7) == "sgData|") {
        if (text.length > 8) {

            // convert  ".]" to "]"  and  ",}" to "}"
            var newtext = text.substring(7);
            newtext = newtext.replace(/,]/g, "]").replace(/,}/g, "}");
            wv.sgData = JSON.parse(newtext);
        } else {
            wv.sgData = {};
        }

    // if it starts with "sgFocus|" store the scene graph focus data
    } else if (text.substring(0,8) == "sgFocus|") {
        if (text.length > 9) {
            wv.focus = JSON.parse(text.substring(8));
        } else {
            wv.focus = [];
        }

//        if (wv.curTool.redraw !== undefined) {
//            wv.curTool.redraw();
//        }

        wv.sceneUpd = 1;

    // if it starts with "getCvals|" build the (global) cval array
    } else if (text.substring(0,9) == "getCvals|") {
        if (wv.server == "serveCAPS") {
            if (text.length > 10) {
                cval = JSON.parse(text.substring(9));
            } else {
                cval = new Array;
            }
            wv.cvalStat = -1;

            rebuildTreeWindow(3);
        }

    // if it starts with "newCval|" do nothing
    } else if (text.substring(0,8) == "newCval|") {

    // if it starts with "setCval|", post message if error detected
    } else if (text.substring(0,8) == "setCval|") {
        if (text.substring(0,15) == "setCval|ERROR::") {
            alert("value(s) not changed because error was detected");

            browserToServer("getCvals|2|");
            wv.cvalStat = 6000;
        }

    // if it starts with "getPmtrs|" build the (global) pmtr array
    } else if (text.substring(0,9) == "getPmtrs|") {
        if (text.length > 10) {
            pmtr = JSON.parse(text.substring(9));
        } else {
            pmtr = new Array;
        }
        wv.pmtrStat = -1;

        rebuildTreeWindow(4);

    // if it starts with "newPmtr|" do nothing
    } else if (text.substring(0,8) == "newPmtr|") {

    // if it starts with "clrVels|" do nothing
    } else if (text.substring(0,8) == "clrVels|") {

    // if it starts with "setPmtr|", post message if error detected
    } else if (text.substring(0,8) == "setPmtr|") {
        if (text.substring(0,15) == "setPmtr|ERROR::") {
            alert("value(s) not changed because error was detected");

            browserToServer("getPmtrs|");
            wv.pmtrStat = 6000;
        }

    // if it starts with "delValu|" do nothing
    } else if (text.substring(0,8) == "delValu|") {

    // if it starts with "setVel|" do nothing
    } else if (text.substring(0,7) == "setVel|") {

    // if it starts with "getBrchs|" build the (global) brch array
    } else if (text.substring(0,9) == "getBrchs|") {
        if (wv.server != "serveCAPS") {
            if (text.length > 10) {
                brch = JSON.parse(text.substring(9));
            } else {
                brch = new Array;
            }
            wv.brchStat = -1;

            rebuildTreeWindow(5);
        }

    // if it starts with "newBrch|" do nothing (except possibly post warning)
    } else if (text.substring(0,8) == "newBrch|") {
        if (text.substring(8,17) == "WARNING::") {
            postMessage(text.substring(8));
        }

    // if it starts with "setBrch|" do nothing (except possibly post warning)
    } else if (text.substring(0,8) == "setBrch|") {
        if (text.substring(8,17) == "WARNING::") {
            postMessage(text.substring(8));
        }

    // if it starts with "delBrch|" do nothing (except possibly post warning)
    } else if (text.substring(0,8) == "delBrch|") {
        if (text.substring(8,17) == "WARNING::") {
            postMessage(text.substring(8));
        }

    // if it starts with "setAttr|" do nothing
    } else if (text.substring(0,8) == "setAttr|") {

    // if it starts with "undo|" get the latest Parameters and
    //                           Branches and update tree (scene)
    } else if (text.substring(0,5) == "undo|") {
        var cmd = text.substring(5,text.length-1);

        wv.sgUpdate = 1;
        wv.brchStat = 0;
        wv.pmtrStat = 0;
        wv.cvalStat = 0;

        postMessage("Undoing '"+cmd+"' ====> Re-build is needed <====");

        activateBuildButton();

    // if it starts with "new|" post message
    } else if (text.substring(0,4) == "new|") {
        postMessage("New (blank) configuration has been loaded\n" +
                    "    Press \"Design Parameters\" (in left window) to add a Design Parameter\n" +
                    "    Press \"Branches\" (in left window) to begin a 3D Object\n" +
                    "    Press \"Branches\" and then \"skbeg\" to begin a 2D Sketch" +
                    "    Press \"File->Open\" to open an existing file\n" +
                    "    Press \"File->Edit<new file>\" to edit a new file");

        caps = new Array();
        pmtr = new Array();
        brch = new Array();

        rebuildTreeWindow(6);

    // if it starts with "save|" do nothing
    } else if (text.substring(0,5) == "save|") {

    // if it starts with "load|", ask for a build
    } else if (text.substring(0,5) == "load|") {

        // if only one file, inform the user that a rebuild is in process
        if (wv.filenames.split("|").length <= 3) {
            var button = document.getElementById("solveButton");
            button["innerHTML"] = "Re-building...";
            button.style.backgroundColor = "#FFFF3F";

            // turn the background of the message window back to original color
            var botm = document.getElementById("brframe");
            botm.style.backgroundColor = "#F7F7F7";

            browserToServer("build|0|");

            if (wv.server == "serveCAPS") {
                browserToServer("getCvals|3|");
                wv.cvalStat = 6000;
            }

            browserToServer("getPmtrs|");
            wv.pmtrStat = 6000;

            if (wv.server != "serveCAPS") {
                browserToServer("getBrchs|2|");
                wv.brchStat = 6000;
            }

            // inactivate buttons until build is done
            changeMode( 0);
            changeMode(-1);

        // otherwise, activate the Re-Build button
        } else {
            changeMode(0);
            activateBuildButton();
            postMessage("====> Re-build is needed <====");
        }

    // if it starts with "build|" reset the build button
    } else if (text.substring(0,6) == "build|") {
        var button = document.getElementById("solveButton");
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = null;

        // turn the background of the message window back to original color
        var botm = document.getElementById("brframe");
        botm.style.backgroundColor = "#F7F7F7";

        var textList = text.split("|");

        // post messages sent from OpenCSM
        if (textList[3].length > 0) {
            postMessage(textList[3]);
        }

        if (textList[1].substring(0,7) == "ERROR::") {
            alert("how did we get here???");
//            postMessage(textList[1]);
//
//            var ibrch    = Number(textList[2]);
//            var nbody    = Number(textList[3]);
//
//            wv.builtTo = ibrch;
//
//            if (nbody > 0) {
//                postMessage("Build complete through ibrch="+ibrch+
//                            ", which generated "+nbody+" Body(s)");
//            } else {
//                alert("No Bodys were produced");
//            }
        } else {
            var ibrch  = Number(textList[1]);
            var nbody  = Number(textList[2]);

            wv.builtTo = ibrch;

            if (ibrch == 0 && brch.length == 0) {
                // post nothing because we started without a file
            } else if (ibrch == brch.length || brch.length == 0) {
                postMessage("\nEntire build complete, which generated "+nbody+
                            " Body(s)");
                changeMode(0);
            } else if (ibrch >= brch.length) {
                postMessage("\nBuild complete through ibrch="+ibrch+
                            ", which generated "+nbody+" Body(s)");
                changeMode(0);
            } else if (ibrch > 0) {
                postMessage("\nPartial build (through "+brch[ibrch-1].name+
                            ") complete, which generated "+nbody+" Body(s)");
                changeMode(0);
            } else {
//                postMessage("ibrch="+ibrch+"   brch.length="+brch.length);
                postMessage("Build failed before first Branch");
                changeMode(0);
            }
        }

    // if it starts with "loadSketch|" initialize the Sketcher
    } else if (text.substring(0,11) == "loadSketch|") {

        var textList = text.split("|");

        if (textList[1].substring(0,7) == "ERROR::") {
            postMessage(textList[1]);
            alert("Sketch cannot be loaded");
        } else {
            // open the Sketcher
            changeMode(8);
            sketcherInitialize();
            sketcherLoad(textList[1], textList[2], textList[3], textList[4]);
        }

    // if it starts with "solveSketch|" update the Sketcher
    } else if (text.substring(0,12) == "solveSketch|") {

        sketcherSolvePost(text);

    // if it starts with "saveSketchBeg|" do nothing
    } else if (text.substring(0,14) == "saveSketchBeg|") {

    // if it starts with "saveSketchMid|" do nothing
    } else if (text.substring(0,14) == "saveSketchMid|") {

    // if it starts with "saveSketchEnd|" do nothing
    } else if (text.substring(0,14) == "saveSketchEnd|") {

    // if it starts with "saveSketch|" either return an error or update the mode
    } else if (text.substring(0,11) == "saveSketch|") {

        if (text.substring(11,14) == "ok|") {
            activateBuildButton();
            changeMode(0);
        } else {
            alert("The sketch could not be saved");
        }

    // if it starts with "makeEBody|", report the error if there is one
    } else if (text.substring(0,10) == "makeEBody|") {
        if (text.substring(10,17) == "ERROR::") {
            alert(text.substring(17));

            var button = document.getElementById("solveButton");
            button["innerHTML"] = "ShowEBodys";
        }

    // if it starts with "getBodyDetails|", post the reposnse
    } else if (text.substring(0,15) == "getBodyDetails|") {
        var textList = text.split("|");

        postMessage(textList[3]);
        postMessage(" ");

    // if it starts with "caps|list|" post a message */
    } else if (text.substring(0,9) == "caps|list") {
        var textList = text.split("|");

        postMessage(textList[2]);
        postMessage(" ");
        
    // if it starts with "caps|phase|" initialize the Phase */
    } else if (text.substring(0,5) == "caps|") {
        var textList = text.split("|");

        /* get the Cvals and Pmtrs */
        browserToServer("getCvals|0|");
        wv.cvalStat = -1;

        browserToServer("getPmtrs|");
        wv.pmtrStat = -1;

        wv.capsPhase = textList[1];

        wv.capsIntent = prompt("Working on Phase "+wv.capsPhase+" of Project \""+wv.capsProj+"\"\n"+
                               "What is the intent of the Phase?");
        
    // if it starts with "timLoad|" pass to curTool or postMessage
    } else if (text.substring(0,8) == "timLoad|") {
        if (wv.overlay !== undefined && wv.overlay.timLoadCB !== undefined) {
            wv.overlay.timLoadCB(text.substring(8));
        } else if (wv.curTool.timLoadCB !== undefined) {
            wv.curTool.timLoadCB(text.substring(8));
        } else {
            postMessage(text);
        }

    // if it starts with "timSave|" pass to curTool or postMessage
    } else if (text.substring(0,8) == "timSave|") {
        if (wv.overlay !== undefined && wv.overlay.timSaveCB !== undefined) {
            wv.overlay.timSaveCB(text.substring(8));
        } else if (wv.curTool.timSaveCB !== undefined) {
            wv.curTool.timSaveCB(text.substring(8));
        } else {
            postMessage(text);
        }

    // if it starts with "timQuit|" pass to curTool or postMessage
    } else if (text.substring(0,8) == "timQuit|") {
//$$$        alert("in wvServerMessage: wv.overlay="+wv.overlay+"   wv.curTool.name="+wv.curTool.name+"   text="+text);
        if (wv.overlay !== undefined && wv.overlay.timQuicBB !== undefined) {
            wv.overlay.timQuitCB(text.substring(8));
        } else if (wv.curTool.timQuitCB !== undefined) {
//$$$            alert("here");
            wv.curTool.timQuitCB(text.substring(8));
        } else {
//$$$            alert("there");
            postMessage(text);
        }

    // if it starts with "timMesg|" pass to curTool or postMessage
    } else if (text.substring(0,8) == "timMesg|") {
        if (wv.overlay !== undefined && wv.overlay.timMesgCB !== undefined) {
            wv.overlay.timMesgCB(text.substring(8));
        } else if (wv.curTool.timMesgCB !== undefined) {
            wv.curTool.timMesgCB(text.substring(8));
        } else {
            postMessage(text);
        }

    // if it starts with "timDraw|" do nothing
    } else if (text.substring(0,8) == "timDraw|") {

    // if it starts with "overlayBeg|" start overlay mode
    } else if (text.substring(0,11) == "overlayBeg|") {

        cmdOverlayBeg(text.substring(11));

    // if it starts with "postMessage|", post the message
    } else if (text.substring(0,12) == "postMessage|") {
        postMessage(text.substring(12));

    // if it starts with "setLims|" do nothing
    } else if (text.substring(0,8) == "setLims|") {

    // if it matches "getFilenames||", tell user how to begin
    } else if (text == "getFilenames||") {
        if (wv.fileintro === false) {
            wv.fileintro = true;

            postMessage("ESP has started without a .csm file\n" +
                        "    Press \"Design Parameters\" (in left window) to add a Design Parameter\n" +
                        "    Press \"Branches\" (in left window) to begin a 3D Object\n" +
                        "    Press \"Branches\" and then \"skbeg\" to begin a 2D Sketch\n" +
                        "    Press \"File->Open\" to open an existing file\n" +
                        "    Press \"File->Edit<new file>\" to edit a new file");
        }

    // if it starts with "getFilenames|" store the results in wv.filenames
    } else if (text.substring(0,13) == "getFilenames|") {
        var oldName = "";
        if (wv.filenames.length != 1) {
            var foo = wv.filenames.split("|");
            oldName = foo[1];
        }

        wv.filenames = text.substring(12);
        var textList = text.split("|");

        // if in CAPS mode, change reference from capsRestart.cpc to capsProj.csm
        if (wv.server == "serveCAPS") {
            var newText = textList[1].replace("/"+wv.capsPhase+"/capsRestart.cpc",
                                              "/"+wv.capsProj+".csm");
            textList[1]  = newText;
            wv.filenames = textList[1];
            for (var i = 2; i < textList.length-1; i++) {
                wv.filenames += "|" + textList[i];
            }
        }

        if (textList[1].length > 0 && textList[1] == oldName) {
            var mesg = "";
        } else if (textList[1].length > 0) {
            var mesg = "\"" + textList[1] + "\" has been loaded";

            for (var i = 2; i < textList.length-1; i++) {
                if (textList[i].length > 0) {
                    mesg += "\n    uses: \"" + textList[i] + "\"";
                }
            }
        } else {
            var mesg = "no .csm file has been loaded";
        }
        postMessage(mesg);

    // if it starts with "getCsmFile|" store the results in wv.curFile and change to mode 7
    } else if (text.substring(0,11) == "getCsmFile|") {
        wv.curFile = text.substring(11);

        // special treatment if entering Gloves
        if (wv.curMode == 9) {

            // load Gloves from wv.curFile
            glovesLoad();

            // if wv.curFile does not have a Gloves section, ask for a component
            if (gloves.comp.length == 0) {
                gloves.cmdSolve();
            }

            // if there are still no Gloves components, quit Gloves
            if (gloves.comp.length == 0) {
                gloves.cmdQuit();

            // if there is at least one component, show it/them
            } else {
                glovesDraw();
                gloves.cmdHome();
            }

            return;
        }

        // remove former editor
        if (wv.codeMirror !== undefined && wv.myRole != 0) {
            wv.codeMirror.toTextArea();
            wv.codeMirror = undefined;
            wv.fileindx   = undefined;
        }

        // fill in the name of the .csm file
        var editorFilename = document.getElementById("editorFilename");
        if (wv.myRole != 0) {
            var filelist = wv.filenames.split("|");
            if (wv.fileindx === undefined && filelist.length == 3) {
                wv.fileindx = 1;
            }
            editorFilename["innerHTML"] = "Contents of: "+filelist[wv.fileindx]+" (read-only)";
        } else if (wv.fileindx >= 0) {
            var filelist = wv.filenames.split("|");
            editorFilename["innerHTML"] = "Contents of: "+filelist[wv.fileindx];
        } else {
            editorFilename["innerHTML"] = "Contents of: &lt new file &gt";
        }

        // fill the textarea with the current .csm file
        var editorTextArea = document.getElementById("editorTextArea");

        editorTextArea.cols  = 84;
        editorTextArea.rows  = 25;
        editorTextArea.value = wv.curFile;

        // unhide so that CodeMirror initialization will work
        var editorForm    = document.getElementById("editorForm");
        editorForm.hidden = false;

        // initialize CodeMirror
        if (wv.codeMirror === undefined) {
            var options = {
                lineNumbers:                     true,
                mode:                      "csm_mode",
                smartIndent:                     true,
                indentUnit:                         3,
                scrollbarStyle:              "native",
                matchBrackets:                   true,
                styleActiveLine:     {nonEmpty: true},
                styleActiveSelected: {nonEmpty: true},
                theme:                       "simple"
            };

            wv.codeMirror = CodeMirror.fromTextArea(editorTextArea, options);

            // choose between editable and readonly
            if (wv.myRole == 0) {
                wv.codeMirror.setOption("readOnly", false);
            } else {
                wv.codeMirror.setOption("readOnly",  true);
            }

            // send changes made in the editor to other collaborative users
            CodeMirror.on(wv.codeMirror, "change",
                          function (instance, obj) {
                              if (wv.myRole == 0 && wv.numUsers > 1) {
                                  if (obj.removed == "") {
                                      if (obj.text.length == 2) {
                                          browserToServer("editor|add|"+obj.from.line+"|"+obj.from.ch+"|<cr>|");
                                      } else {
                                          browserToServer("editor|add|"+obj.from.line+"|"+obj.from.ch+"|"+obj.text+"|");
                                      }
                                  } else {
                                      browserToServer("editor|del|"+obj.from.line+"|"+obj.from.ch+"|"+obj.to.line+"|"+obj.to.ch+"|");
                                  }
                              }
                          });

            var topRite = document.getElementById("trframe");
            var height  = Number(topRite.style.height.replace("px",""));
            wv.codeMirror.setSize(null, height-65);  // good value for nominal browser fonts
        }

        // clear any previous undo history and clipboard
        wv.codeMirror.clearHistory();
        wv.clipboard = "";

        // make wv.linenum the active line

        if (wv.linenum != 0) {
            wv.codeMirror.focus();
            wv.codeMirror.setCursor({line: wv.linenum-1, ch: 0});
            wv.linenum = 0;

        // if we previously editted this file, go back to line number
        } else if (wv.fileindx >= 0) {
            var filelist = wv.filenames.split("|");
            if (filelist[wv.fileindx] == wv.lastfile && wv.lastline >= 0) {
                wv.codeMirror.focus();
                wv.codeMirror.setCursor({line: wv.lastline-1, ch: 0});
            }
        }

        wv.lastfile = "";
        wv.lastline = -1;

        // post the editorForm
        changeMode(7);

    // if it starts with "setCsmFileBeg|" do nothing
    } else if (text.substring(0,14) == "setCsmFileBeg|") {

    // if it starts with "setCsmFileMid|" do nothing
    } else if (text.substring(0,14) == "setCsmFileMid|") {

    // if it starts with "setCsmFileEnd|" do nothing
    } else if (text.substring(0,14) == "setCsmFileEnd|") {

    // if it starts with "setWvKey|" turn key or logo on
    } else if (text.substring(0,9) == "setWvKey|") {
        if (wv.curMode == 0 || wv.curMode == 12) {
            if (text.substring(9,11) == "on") {
                document.getElementById("WVkey"  ).hidden = false;
                document.getElementById("ESPlogo").hidden = true;
            } else {
                document.getElementById("WVkey"  ).hidden = true;
                document.getElementById("ESPlogo").hidden = false;
            }
        }

    // if it starts with "saveView|" do nothing
    } else if (text.substring(0,9) == "saveView|") {

    // if it starts with "readView|" load view matrix and update display
    } else if (text.substring(0,9) == "readView|") {
        var readViewList = text.split("|");
        var entries      = readViewList[2].split(",");
        if (entries.length == 16) {
            var matrix  = Array(16);
            for (var i = 0; i < 16; i++) {
                matrix[i] = parseFloat(entries[i]);
            }
            wv.mvMatrix.makeIdentity();
            wv.uiMatrix.load(matrix);
            wv.scale    = parseFloat(readViewList[1]);
            wv.sceneUpd = 1;
        } else {
            postMessage("File has wrong number of entries");
        }

    // if it starts with "editor|" and you do not have the ball, modify the text area
    } else if (text.substring(0,7) == "editor|") {
        if (wv.myRole != 0) {
            var textList = text.split("|");

            if        (textList[1] == "add" && textList[4] == "<cr>") {
                wv.codeMirror.getDoc().replaceRange("\n",
                                                    {line: Number(textList[2]), ch: Number(textList[3])});
            } else if (textList[1] == "add") {
                wv.codeMirror.getDoc().replaceRange(textList[4],
                                                    {line: Number(textList[2]), ch: Number(textList[3])});
            } else if (textList[1] == "del") {
                wv.codeMirror.getDoc().replaceRange("",
                                                    {line: Number(textList[2]), ch: Number(textList[3])},
                                                    {line: Number(textList[4]), ch: Number(textList[5])});
            }
        }

    // if it starts with "xform|" and we are sync'd, load matrix and update display
    } else if (text.substring(0,6) == "xform|") {
        if (wv.myRole == 1) {
            var textList     = text.split("|");
            var masterWidth  = Number(textList[1]);
            var masterHeight = Number(textList[2]);
            var entries      = textList[4].split(",");
            if (entries.length == 16) {
                var matrix  = Array(16);
                for (var i = 0; i < 16; i++) {
                    matrix[i] = parseFloat(entries[i]);
                }
                wv.mvMatrix.makeIdentity();
                wv.uiMatrix.load(matrix);
                wv.scale    = parseFloat(textList[3]);
                wv.sceneUpd = 1;

//                // account for different aspect ratios in browsers
//                var xscale = wv.width  / masterWidth;
//                var yscale = wv.height / masterHeight;
//
//                if (xscale < yscale) {
//                    wv.uiMatrix.scale(xscale, xscale, xscale);
//                } else {
//                    wv.uiMatrix.scale(yscale, yscale, yscale);
//                }
            }
            wv.lastXform = text;
        } else if (wv.myRole == 2) {
            wv.lastXform = text;
        }

    // if it starts with "lastPoint|" and we are sync'd, update or clear the lastPoint
    } else if (text.substring(0,10) == "lastPoint|") {
        if (wv.myRole == 1) {
            var textList = text.split("|");
            if (textList[1] == "off") {
                wv.lastPoint = undefined;
            } else {
                wv.lastPoint = [Number(textList[1]), Number(textList[2]), Number(textList[3])];
            }
            wv.sceneUpd = 1;
        }

    // if it starts with "resetMode|" and we do not have the ball,
    // reset the mode to 0
    } else if (text.substring(0,10) == "resetMode|") {
        changeMode(0);

    // if it starts with "toggle|" and we are syncd, update the graphic attribute
    } else if (text.substring(0,7) == "toggle|") {
        if (wv.myRole == 1) {
            textList = text.split("|");

            var inode = Number(textList[1]);
            var icol  = Number(textList[2]);
            var state =        textList[3];

            if        (inode == -1) {
                myTree.delay = 1;
            } else if (inode == -2) {
                myTree.delay = 0;
            } else {
                myTree.prop(inode, icol, state);
            }

            myTree.update();
        }

    // if it is "sendState|" and I have the ball, broadcast current
    // settings to everyone
    } else if (text.substring(0,10) == "sendState|") {
        if (wv.myRole == 0) {

            // tell synced browser to delay updating their trees
            browserToServer("toggle|-1|0|0|");

            for (var inode = 0; inode < myTree.parent.length; inode++) {

                // skip if one of the headers
                if (myTree.prop3[inode] === undefined) {
                    continue;
                } else if (myTree.name[inode] == "Display") {
                    continue;
                } else if (myTree.name[inode] == "\u00a0\u00a0\u00a0\u00a0Faces") {
                    continue;
                } else if (myTree.name[inode] == "\u00a0\u00a0\u00a0\u00a0Edges") {
                    continue;
                } else if (myTree.name[inode] == "\u00a0\u00a0\u00a0\u00a0Nodes") {
                    continue;
                } else if (myTree.name[inode] == "\u00a0\u00a0\u00a0\u00a0Csystems") {
                    continue;
                }

                // process this node
                if (myTree.prop1[inode] == "Viz") {
                    var thisNode = document.getElementById("node"+inode+"col3");
                    if (thisNode.getAttribute("class") == "fakelinkon") {
                        changeProp(inode, 1, "on");
                    } else {
                        changeProp(inode, 1, "off");
                    }
                }
                if (myTree.prop2[inode] == "Grd") {
                    var thisNode = document.getElementById("node"+inode+"col4");
                    if (thisNode.getAttribute("class") == "fakelinkon") {
                        changeProp(inode, 2, "on");
                    } else {
                        changeProp(inode, 2, "off");
                    }
                }
                if (myTree.prop3[inode] == "Trn") {
                    var thisNode = document.getElementById("node"+inode+"col5");
                    if (thisNode.getAttribute("class") == "fakelinkon") {
                        changeProp(inode, 3, "on");
                    } else {
                        changeProp(inode, 3, "off");
                    }
                } else if (myTree.prop3[inode] == "Ori") {
                    var thisNode = document.getElementById("node"+inode+"col5");
                    if (thisNode.getAttribute("class") == "fakelinkon") {
                        changeProp(inode, 3, "on");
                    } else {
                        changeProp(inode, 3, "off");
                    }
                }
            }

            // tell synced browsers to update their trees now
            browserToServer("toggle|-2|0|0|");

            // send current transform
            browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
        }

    // if it starts with "nextStep|" either switch button legend back
    //  to StepThru or post message about what is showing
    } else if (text.substring(0,9) == "nextStep|") {
        var nextStepList = text.split("|");
        if (nextStepList[1] == 0) {
            wv.curStep = 0;
            var button = document.getElementById("stepThruBtn");
            button["innerHTML"] = "StepThru";

            postMessage("Finished with StepThru mode");
        } else {
            wv.curStep = Number(nextStepList[1]);
            postMessage("Showing \""+nextStepList[2]+"\" generated by \""+
                       nextStepList[3]+"\" in StepThru mode");
        }

    // if it starts with "ERROR:: there is nothing to undo" post the error
    } else if (text.substring(0,32) == "ERROR:: there is nothing to undo") {
        alert("There is nothing to undo");

    // if it starts with "ERROR:: save(" post the error
    } else if (text.substring(0,13) == "ERROR:: save(") {
        postMessage(text);
        alert("File not written");

    // if it starts with "ERROR::" post the error
    } else if (text.substring(0,7) == "ERROR::") {
        var textList = text.split("|");

        // post messages sent from OpenCSM
        if (textList[1].length > 0) {
            postMessage(textList[1]);
        }

        postMessage(textList[0]);

        // turn the background of the message window back to light-yellow
        var botm = document.getElementById("brframe");
        botm.style.backgroundColor = "#FFFF99";

        var button = document.getElementById("solveButton");

        button["innerHTML"] = "Fix before re-build";
        button.style.backgroundColor = "#FF3F3F";

        changeMode(0);

    // if it starts with "message|" post the message
    } else if (text.substring(0,8) == "message|") {
        postMessage(text.substring(8,text.length-1));

    // default is to post the message
    } else {
        postMessage(text);
    }
};


//
// callback when the server goes down either on abort or normal closure (called by wv-socket.js)
//
var wvServerDown = function () {
    // deactivate the buttons
    wv.curMode = -1;

    postMessage("The server has terminated or network connection has been lost.\n"
                +"ESP will be working off-line");

    // turn the background of the message window pink
    var botm = document.getElementById("brframe");
    botm.style.backgroundColor = "#FFAFAF";
};


//
// callback used to put axes on the canvas (called by ESP.html)
//
var wvUpdateCanvas = function (gl) {

    // draw the lastPoint if it has been defined
    drawLastPoint(gl);

    // draw any wireframes that might exist
    drawWireframes(gl);

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
    var alen =  0.25;     // length of axes

    // set up coordinates for axes
    mVal = mv.getAsArray();

    var verts = new Float32Array(66);
    verts[ 0] = x;
    verts[ 1] = y;
    verts[ 2] = z;
    verts[ 3] = x + alen*(    mVal[ 0]             );
    verts[ 4] = y + alen*(    mVal[ 1]             );
    verts[ 5] = z + alen*(    mVal[ 2]             );
    verts[ 6] = x + alen*(1.1*mVal[ 0]+0.1*mVal[ 4]);
    verts[ 7] = y + alen*(1.1*mVal[ 1]+0.1*mVal[ 5]);
    verts[ 8] = z + alen*(1.1*mVal[ 2]+0.1*mVal[ 6]);
    verts[ 9] = x + alen*(1.3*mVal[ 0]-0.1*mVal[ 4]);
    verts[10] = y + alen*(1.3*mVal[ 1]-0.1*mVal[ 5]);
    verts[11] = z + alen*(1.3*mVal[ 2]-0.1*mVal[ 6]);
    verts[12] = x + alen*(1.1*mVal[ 0]-0.1*mVal[ 4]);
    verts[13] = y + alen*(1.1*mVal[ 1]-0.1*mVal[ 5]);
    verts[14] = z + alen*(1.1*mVal[ 2]-0.1*mVal[ 6]);
    verts[15] = x + alen*(1.3*mVal[ 0]+0.1*mVal[ 4]);
    verts[16] = y + alen*(1.3*mVal[ 1]+0.1*mVal[ 5]);
    verts[17] = z + alen*(1.3*mVal[ 2]+0.1*mVal[ 6]);

    verts[18] = x;
    verts[19] = y;
    verts[20] = z;
    verts[21] = x + alen*(    mVal[ 4]             );
    verts[22] = y + alen*(    mVal[ 5]             );
    verts[23] = z + alen*(    mVal[ 6]             );
    verts[24] = x + alen*(1.1*mVal[ 4]+0.1*mVal[ 8]);
    verts[25] = y + alen*(1.1*mVal[ 5]+0.1*mVal[ 9]);
    verts[26] = z + alen*(1.1*mVal[ 6]+0.1*mVal[10]);
    verts[27] = x + alen*(1.2*mVal[ 4]             );
    verts[28] = y + alen*(1.2*mVal[ 5]             );
    verts[29] = z + alen*(1.2*mVal[ 6]             );
    verts[30] = x + alen*(1.3*mVal[ 4]+0.1*mVal[ 8]);
    verts[31] = y + alen*(1.3*mVal[ 5]+0.1*mVal[ 9]);
    verts[32] = z + alen*(1.3*mVal[ 6]+0.1*mVal[10]);
    verts[33] = x + alen*(1.2*mVal[ 4]             );
    verts[34] = y + alen*(1.2*mVal[ 5]             );
    verts[35] = z + alen*(1.2*mVal[ 6]             );
    verts[36] = x + alen*(1.2*mVal[ 4]             );
    verts[37] = y + alen*(1.2*mVal[ 5]             );
    verts[38] = z + alen*(1.2*mVal[ 6]             );
    verts[39] = x + alen*(1.2*mVal[ 4]-0.1*mVal[ 8]);
    verts[40] = y + alen*(1.2*mVal[ 5]-0.1*mVal[ 9]);
    verts[41] = z + alen*(1.2*mVal[ 6]-0.1*mVal[10]);

    verts[42] = x;
    verts[43] = y;
    verts[44] = z;
    verts[45] = x + alen*(    mVal[ 8]             );
    verts[46] = y + alen*(    mVal[ 9]             );
    verts[47] = z + alen*(    mVal[10]             );
    verts[48] = x + alen*(1.1*mVal[ 8]+0.1*mVal[ 0]);
    verts[49] = y + alen*(1.1*mVal[ 9]+0.1*mVal[ 1]);
    verts[50] = z + alen*(1.1*mVal[10]+0.1*mVal[ 2]);
    verts[51] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
    verts[52] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
    verts[53] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
    verts[54] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
    verts[55] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
    verts[56] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
    verts[57] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
    verts[58] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
    verts[59] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
    verts[60] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
    verts[61] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
    verts[62] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
    verts[63] = x + alen*(1.3*mVal[ 8]-0.1*mVal[ 0]);
    verts[64] = y + alen*(1.3*mVal[ 9]-0.1*mVal[ 1]);
    verts[65] = z + alen*(1.3*mVal[10]-0.1*mVal[ 2]);

    // set up colors for the axes
    var colors = new Uint8Array(66);
    colors[ 0] = 255;   colors[ 1] =   0;   colors[ 2] =   0;
    colors[ 3] = 255;   colors[ 4] =   0;   colors[ 5] =   0;
    colors[ 6] = 255;   colors[ 7] =   0;   colors[ 8] =   0;
    colors[ 9] = 255;   colors[10] =   0;   colors[11] =   0;
    colors[12] = 255;   colors[13] =   0;   colors[14] =   0;
    colors[15] = 255;   colors[16] =   0;   colors[17] =   0;
    colors[18] =   0;   colors[19] = 255;   colors[20] =   0;
    colors[21] =   0;   colors[22] = 255;   colors[23] =   0;
    colors[24] =   0;   colors[25] = 255;   colors[26] =   0;
    colors[27] =   0;   colors[28] = 255;   colors[29] =   0;
    colors[30] =   0;   colors[31] = 255;   colors[32] =   0;
    colors[33] =   0;   colors[34] = 255;   colors[35] =   0;
    colors[36] =   0;   colors[37] = 255;   colors[38] =   0;
    colors[39] =   0;   colors[40] = 255;   colors[41] =   0;
    colors[42] =   0;   colors[43] =   0;   colors[44] = 255;
    colors[45] =   0;   colors[46] =   0;   colors[47] = 255;
    colors[48] =   0;   colors[49] =   0;   colors[50] = 255;
    colors[51] =   0;   colors[52] =   0;   colors[53] = 255;
    colors[54] =   0;   colors[55] =   0;   colors[56] = 255;
    colors[57] =   0;   colors[58] =   0;   colors[59] = 255;
    colors[60] =   0;   colors[61] =   0;   colors[62] = 255;
    colors[63] =   0;   colors[64] =   0;   colors[65] = 255;

    // draw the axes
    if (gl.lineWidth) {
        gl.lineWidth(2);
    }
    gl.disableVertexAttribArray(2);
    gl.uniform1f(wv.u_wLightLoc, 0.0);

    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, verts, gl.STATIC_DRAW);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(0);

    var cbuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
    gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
    gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.enableVertexAttribArray(1);

    gl.drawArrays(gl.LINES, 0, 22);
    gl.deleteBuffer(buffer);
    gl.deleteBuffer(cbuf);
    gl.uniform1f(wv.u_wLightLoc, 1.0);
};


//
// callback when the scene has been updated (called by wv-socket.js)
//
wv.sceneUpdated = function () {
    // alert ("in wv.sceneUpdated()");

    if (wv.curTool.sceneUpdated !== undefined) {
        wv.curTool.sceneUpdated();
    }
};

////////////////////////////////////////////////////////////////////////////////


//
// callback when "fileButton" is pressed (called by ESP.html)
//
var cmdFile = function () {
    // alert("in cmdFile()");

    if (checkIfWithBall() === false) {
        return;
    }

    // toggle between hiding and showing the File menu contents
    document.getElementById("myFileMenu").classList.toggle("showFileMenu");

    // if myToolMenu is currently posted, delete it now
    document.getElementById("myToolMenu").classList.remove("showToolMenu");

    // if myCollabMenu is currently posted, delete it now
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    /* remove previous menu entries */
    var menu = document.getElementById("myFileMenu");
    while (menu.firstChild !== null) {
        menu.removeChild(menu.firstChild);
    }

    /* add entries for New, Open, Export, and Edit for all files */
    var button;

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Create new model";
    button.value   = "New";
    button.onclick = cmdFileNew;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Open .csm file";
    button.value   = "Open";
    button.onclick = cmdFileOpen;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Export .csm file from feature tree";
    button.value   = "Export FeatureTree";
    button.onclick = cmdFileExport;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    if (wv.filenames.length > 0) {
        var filelist = wv.filenames.split("|");
        for (var ielem = 0; ielem < filelist.length; ielem++) {
            if (filelist[ielem].length > 0) {
                button = document.createElement("input");
                button.type     = "button";
                button.title    = "Edit the file";
                button.value    = "Edit: "+filelist[ielem];
                button.fileindx = ielem;
                button.onclick  = function (ielem) {cmdFileEdit(ielem);};

                menu.appendChild(document.createElement("br"));
                menu.appendChild(button);
            }
        }
    }

    button = document.createElement("input");
    button.type     = "button";
    button.title    = "Edit new file";
    button.value    = "Edit: <new file>";
    button.fileindx = -1;
    button.onclick  = function (ielem) {cmdFileEdit(ielem);};
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);
};


//
// callback when "File->New" is pressed
//
var cmdFileNew = function () {
    // alert("in cmdFileNew()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.nchanges > 0) {
        if (confirm(wv.nchanges+" change(s) will be lost.  Continue?") !== true) {
            return;
        }
    } else {
        if (confirm("Do you really want to remove the current model?") !== true) {
            return;
        }
    }

    if (wv.curMode == 0) {
        browserToServer("new|");

        wv.filenames = "";
        wv.nchanges  = 0;

        // remove from editor
        if (wv.codeMirror !== undefined) {
            wv.codeMirror.toTextArea();
            wv.codeMirror = undefined;
        }

        cval   = new Array();
        pmtr   = new Array();
        brch   = new Array();
        sgData = {};
    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
};


//
// callback when "File->Open" is pressed
//
var cmdFileOpen = function () {
    // alert("in cmdFileOpen()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.nchanges > 0) {
        if (confirm(wv.nchanges+" change(s) will be lost.  Continue?") !== true) {
            return;
        }
    }

    if (wv.curMode == 0) {
        var filelist = wv.filenames.split("|");
        var filename = prompt("Enter filename to open:", filelist[1]);
        if (filename !== null) {
            if (filename.length == 0) {
                alert("empty filename given");
                return;
            } else if (filename.search(/\.csm$/) > 0 ||
                       filename.search(/\.cpc$/) > 0   ) {
                // well-formed filename
            } else if (filename.search(/\.udc$/) > 0   ) {
                alert("cannot open .udc with File->Open");
                return;
            } else {
                // add .csm extension
                filename += ".csm";
            }

            postMessage("Opening \""+filename+"\" ...");

            browserToServer("open|"+filename+"|");

            wv.filenames = filename;
            wv.nchanges  = 0;

            // remove former editor
            if (wv.codeMirror !== undefined) {
                wv.codeMirror.toTextArea();
                wv.codeMirror = undefined;
            }

            if (wv.server == "serveCAPS") {
                browserToServer("getCvals|4|");
                wv.cvalStat = 6000;
            }

            browserToServer("getPmtrs|");
            wv.pmtrStat = 6000;

            if (wv.server != "serveCAPS") {
                browserToServer("getBrchs|3|");
                wv.brchStat = 6000;
            }

            var button = document.getElementById("solveButton");
            button["innerHTML"] = "Re-building...";
            button.style.backgroundColor = "#FFFF3F";

            //inactivate buttons until build is done
            changeMode(-1);
        } else {
            postMessage("NOT opening since no filename specified");
        }

    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
};


//
// callback when "File->ExportFeatureTree" is pressed
//
var cmdFileExport = function () {
    // alert("in cmdFileExport()");

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.curMode == 0) {
        var filelist = wv.filenames.split("|");
        var filename = prompt("Enter filename to write:", filelist[0]);
        if (filename !== null) {
            if (filename.search(/\.csm$/) > 0 ||
                filename.search(/\.cpc$/) > 0 ||
                filename.search(/\.udc$/) > 0   ) {
                // well-formed filename
            } else {
                // add .csm extension
                filename += ".csm";
            }

            // warn that formatting will be lost
            if (confirm("This will export a file by reading the"
                        + " FeatureTree and DesignParameters.\n  As"
                        + " a result, you will lose all your"
                        + " formatting and comments.  Continue?") !== true) {
                return;
            }

            postMessage("Saving model to '"+filename+"'");
            browserToServer("save|"+filename+"|");

            wv.filenames = filename;
            wv.nchanges  = 0;
        } else {
            postMessage("NOT saving since no filename specified");
        }

    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
};


//
// callback when "editButton" is pressed
//
var cmdFileEdit = function (e, indx) {
    // alert("in cmdFileEdit(e="+e+", indx="+indx+")");

    var index;
    if (e !== null) {
        index = e.srcElement.fileindx;
    } else {
        index = indx;
    }

    // close the File menu
    var menu = document.getElementsByClassName("fileMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showFileMenu")) {
            menu[i].classList.remove(  "showFileMenu");
        }
    }

    if (wv.curMode != 0 && e !== null) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.curMode != 0 && indx != wv.fileIndx) {
        var newFile = wv.codeMirror.getDoc().getValue();

        if (wv.curFile == newFile && wv.nchanges == 0) {
            editorCancel();
        } else {
            alert("Current file has been changed, so Cancel current edit first.");
            return;
        }
    }

    // if all you have is a .udc, do not allow editing
    var filenames = wv.filenames.split("|");
    if (filenames[0].includes(".udc") && index >= 0) {
        alert("\""+filenames[0]+"\" cannot be edited because there is no .csm file\n"
              +"use File->New or File->Edit: <new file> to create a .csm file first");
        return;
    }

    // if there have been any changes, tell user to save first
    if (wv.nchanges > 0) {
        if (confirm(wv.nchanges+" change(s) will be lost.  Continue?") !== true) {
            return;
        }
    }

    // create a new file if index < 0
    if (index < 0) {
        browserToServer("getCsmFile|");
        wv.fileindx = -1;
        return;
    }

    // get the latest csm file
    try {
        var filelist = wv.filenames.split("|");
        browserToServer("getCsmFile|"+filelist[index]+"|");
        wv.fileindx = index;
    } catch (e) {
        // could not send message
    }
};


//
// callback when "OK" button is pressed in editorForm (called by ESP.html)
//
var editorOk = function () {
    // alert("in editorOk()");

    if (checkIfWithBall() === false) {
        return;
    }

    // remember what file we are editting and which line we were on
    var myLastline = wv.codeMirror.getCursor()["line"] + 1;

    if (wv.fileindx >= 0) {
        var filelist = wv.filenames.split("|");
        wv.lastfile = filelist[wv.fileindx];
        wv.lastline = myLastline;
    } else {
        wv.lastfile = "";
        wv.lastline = -1;
    }

    // tell user if no changes were made
    var newFile = wv.codeMirror.getDoc().getValue();

    if (wv.curFile == newFile && wv.nchanges == 0 &&
        confirm("No changes were made.  Reload file anyway?") !== true) {

        wv.curFile = "";

        // remove from editor
        if (wv.codeMirror !== undefined) {
            wv.codeMirror.toTextArea();
            wv.codeMirror = undefined;
        }

        changeMode(0);
        return;

    // get filename
    } else if (wv.fileindx < 0 || wv.filenames.length == 0) {
        wv.filenames = prompt("Enter filename");
        if (wv.filenames !== null) {
            if (wv.filenames.search(/\.csm$/) <= 0 &&
                wv.filenames.search(/\.udc$/) <= 0   ) {
                // add .csm extension if no extension is given
                wv.filenames += ".csm";
            }
        } else {
            alert("NOT saving since no filename specified");
            return;
        }

        wv.lastfile = wv.filenames;
        wv.lastline = myLastline;

//$$$    // check for overwrite
//$$$    } else if (confirm("This may overwrite an existing file. " +
//$$$                       "Continue?") !== true) {
//$$$        return;
    }

    // because of an apparent limit on the size of text
    //    messages that can be sent from the browser to the
    //    server, we need to send the new file back in
    //    pieces and then reassemble on the server
    var maxMessageSize = 800;

    var ichar    = 0;
    var part     = newFile.substring(ichar, ichar+maxMessageSize);

    if (wv.fileindx < 0) {
        var myFilename = wv.filenames;
    } else {
        var filelist = wv.filenames.split("|");
        var myFilename = filelist[wv.fileindx];
    }

    browserToServer("setCsmFileBeg|"+myFilename+"|"+part);
    ichar += maxMessageSize;

    while (ichar < newFile.length) {
        part = newFile.substring(ichar, ichar+maxMessageSize);
        browserToServer("setCsmFileMid|"+part);
        ichar += maxMessageSize;
    }

    browserToServer("setCsmFileEnd|");

    // remember the edited .csm file
    wv.curFile = newFile;

    postMessage("'"+myFilename+"' file has been changed.");
    wv.fileindx = -1;

    // get an updated version of the Parameters and Branches
    wv.cvalStat = 0;
    wv.pmtrStat = 0;
    wv.brchStat = 0;

    // remove the contents of the file from memory
    wv.curFile = "";

    // remove from editor
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.toTextArea();
        wv.codeMirror = undefined;
    }

    // reset the number of changes
    wv.nchanges = 0;

    // if only one file, inform the user that a rebuild is in process
    if (wv.filenames.split("|").length <= 3) {
        var button = document.getElementById("solveButton");
        button["innerHTML"] = "Re-building...";
        button.style.backgroundColor = "#FFFF3F";

        // turn the background of the message window back to original color
        var botm = document.getElementById("brframe");
        botm.style.backgroundColor = "#F7F7F7";

        // inactivate buttons until build is done
        changeMode( 0);
        changeMode(-1);

    // otherwise, activate the Re-Build button
    } else {
        changeMode(0);
        activateBuildButton();
        postMessage("====> Re-build is needed <====");
    }
}

var editorOk_old = function () {
    // alert("in editorOk()");

    // tell user if no changes were made
    var newFile = wv.codeMirror.getDoc().getValue();

    if (wv.curFile == newFile && wv.nchanges == 0) {
        alert("No changes were made");
        changeMode(0);
        return;

    // get filename
    } else if (wv.fileindx < 0 || wv.filenames.length == 0) {
        wv.filenames = prompt("Enter filename");
        if (wv.filenames !== null) {
            if (wv.filenames.search(/\.csm$/) <= 0 &&
                wv.filenames.search(/\.udc$/) <= 0   ) {
                // add .csm extension if no extension is given
                wv.filenames += ".csm";
            }
        } else {
            alert("NOT saving since no filename specified");
            return;
        }

    // check for overwrite
    } else if (confirm("This may overwrite an existing file. " +
                       "Continue?") !== true) {
        return;
    }

    // because of an apparent limit on the size of text
    //    messages that can be sent from the browser to the
    //    server, we need to send the new file back in
    //    pieces and then reassemble on the server
    var maxMessageSize = 800;

    var ichar    = 0;
    var part     = newFile.substring(ichar, ichar+maxMessageSize);

    if (wv.fileindx < 0) {
        var myFilename = wv.filenames;
    } else {
        var filelist = wv.filenames.split("|");
        var myFilename = filelist[wv.fileindx];
    }

    browserToServer("setCsmFileBeg|"+myFilename+"|"+part);
    ichar += maxMessageSize;

    while (ichar < newFile.length) {
        part = newFile.substring(ichar, ichar+maxMessageSize);
        browserToServer("setCsmFileMid|"+part);
        ichar += maxMessageSize;
    }

    browserToServer("setCsmFileEnd|");

    // remember the edited .csm file
    wv.curFile = newFile;

    postMessage("'"+myFilename+"' file has been changed.");
    wv.fileindx = -1;

    // get an updated version of the Parameters and Branches
    wv.cvalStat = 0;
    wv.pmtrStat = 0;
    wv.brchStat = 0;

    // remove the contents of the file from memory
    wv.curFile = "";

    // remove from editor
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.toTextArea();
        wv.codeMirror = undefined;
    }

    // reset the number of changes
    wv.nchanges = 0;

    // inform the user that a rebuild is in process
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode( 0);
    changeMode(-1);
};


//
// callback when "Cancel" is pressed in editorForm (called by ESP.html)
//
var editorCancel = function () {
    // alert("in editorCancel()");

    if (checkIfFree() === false) {
        alert("Note: to see the file again, the user with the ball will need to relaunch the editor");
        return;
    }

    // remember what file we are editting and which line we were on
    if (wv.fileindx >= 0) {
        var filelist = wv.filenames.split("|");
        wv.lastfile = filelist[wv.fileindx];
        wv.lastline = wv.codeMirror.getCursor()["line"] + 1;
    } else {
        wv.lastfile = "";
        wv.lastline = -1;
    }

    // remove the contents of the file from memory
    wv.curFile = "";

    // remove from editor
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.toTextArea();
        wv.codeMirror = undefined;
    }

    if (wv.myRole == 0) {
        browserToServer("resetMode|");
    }

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "capsButton" is pressed (called by ESP.html)
//
var cmdCaps = function () {
    // alert("in cmdCaps()");

    if (checkIfWithBall() === false) {
        return;
    }

    // toggle between hiding and showing the Caps menu contents
    document.getElementById("myCapsMenu").classList.toggle("showCapsMenu");

    // if myToolMenu is currently posted, delete it now
    document.getElementById("myToolMenu").classList.remove("showToolMenu");

    // if myCollabMenu is currently posted, delete it now
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    /* remove previous menu entries */
    var menu = document.getElementById("myCapsMenu");
    while (menu.firstChild !== null) {
        menu.removeChild(menu.firstChild);
    }

    /* add entries */
    var button;

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Update intent";
    button.value   = "Update intent";
    button.onclick = cmdCapsUpdateIntent;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Commit Phase";
    button.value   = "Commit Phase";
    button.onclick = cmdCapsCommitPhase;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Revert Phase";
    button.value   = "Revert Phase";
    button.onclick = cmdCapsRevertPhase;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List Phases";
    button.value   = "List Phases";
    button.onclick = cmdCapsListPhases;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List Analyses";
    button.value   = "List Analyses";
    button.onclick = cmdCapsListAnalyses;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List Bounds";
    button.value   = "List Bounds";
    button.onclick = cmdCapsListBounds;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Attach file";
    button.value   = "Attach file";
    button.onclick = cmdCapsAttachFile;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "View files";
    button.value   = "View files";
    button.onclick = cmdCapsViewFiles;
    menu.appendChild(document.createElement("br"));
    menu.appendChild(button);
};


//
// callback when "CAPS->UpdateIntent" is pressed
//
var cmdCapsUpdateIntent = function () {
    // alert("in cmdCapsUpdateIntent()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }

    var intent = prompt("Working on Phase "+wv.capsPhase+" of Project \""+wv.capsProj+"\"\n"+
                        "What is the intent of the Phase?", wv.capsIntent);
    if (intent !== null) {
        if (intent.length > 0) {
            wv.capsIntent = intent;
        }
    }
};


//
// callback when "CAPS->CommitPhase" is pressed
//
var cmdCapsCommitPhase = function () {
    // alert("in cmdCapsCommitPhase()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }

    // close the current Phase and open the next Phase
    browserToServer("caps|commit|");
};


//
// callback when "CAPS->RevertPhase" is pressed
//
var cmdCapsRevertPhase = function () {
    // alert("in cmdCapsRevertPhase()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }

    // close the current Phase (with remove) and reopen
    browserToServer("caps|revert|");
};


//
// callback when "CAPS->ListPhases" is pressed
//
var cmdCapsListPhases = function () {
    // alert("in cmdCapsListPhases()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }

    // get the view
    browserToServer("caps|listPhases|");
};


//
// callback when "CAPS->ListAnalyses" is pressed
//
var cmdCapsListAnalyses = function () {
    // alert("in cmdCapsListAnalyses()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }

    // get the view
    browserToServer("caps|listAnalyses|");
};


//
// callback when "CAPS->ListBounds" is pressed
//
var cmdCapsListBounds = function () {
    // alert("in cmdCapsListBounds()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }

    // get the view
    browserToServer("caps|listBounds|");
};


//
// callback when "CAPS->AttachFile" is pressed
//
var cmdCapsAttachFile = function () {
    alert("in cmdCapsAttachFile()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }
};


//
// callback when "CAPS->ViewFiles" is pressed
//
var cmdCapsViewFiles = function () {
    alert("in cmdCapsViewFiles()");

    // close the Caps menu
    var menu = document.getElementsByClassName("capsMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCapsMenu")) {
            menu[i].classList.remove(  "showCapsMenu");
        }
    }
};


//
// callback when "toolButton" is pressed (called by ESP.html)
//
var cmdTool = function () {
    // alert("in cmdTool()");

    if (checkIfWithBall() === false) {
        return;
    }

    // toggle between hiding and showing the File menu contents
    document.getElementById("myToolMenu").classList.toggle("showToolMenu");

    // if myFileMenu or myCapsMenu is currently posted, delete it now
    document.getElementById("myFileMenu").classList.remove("showFileMenu");
    document.getElementById("myCapsMenu").classList.remove("showCapsMenu");

    // if myCollabMenu is currently posted, delete it now
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    /* remove previous menu entries */
    var menu = document.getElementById("myToolMenu");
    while (menu.firstChild !== null) {
        menu.removeChild(menu.firstChild);
    }

    /* add entries based upon the server */
    var button;

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Launch Sketcher";
    button.value   = "Sketch";
    button.onclick = sketch.launch;
    menu.appendChild(button);

    if (wv.server == "serveESP" || wv.server == "serveCAPS") {
        button = document.createElement("input");
        button.type    = "button";
        button.title   = "Launch ErepEd";
        button.value   = "ErepEd";
        button.onclick = ereped.launch;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "Launch Gloves";
        button.value   = "Gloves";
        button.onclick = gloves.launch;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "Launch Plugs";
        button.value   = "Plugs";
        button.onclick = plugs.launch;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "Launch Pyscript";
        button.value   = "Pyscript";
        button.onclick = pyscript.launch;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "Launch Mitten";
        button.value   = "Mitten";
        button.onclick = mitten.launch;
        menu.appendChild(button);
    }
};


//
// callback when "doneButton" is pressed (called by ESP.html)
//
var cmdDone = function () {
    // alert("in cmdDone()");

    // toggle between hiding and showing the File menu contents
    document.getElementById("myDoneMenu").classList.toggle("showDoneMenu");
};


//
// callback when "stepThruBtn" is pressed (called by ESP.html)
//
var cmdStepThru = function (direction) {
    // alert("in cmdStepThru(direction="+direction+")");

    if (checkIfWithBall() === false) {
        return;
    }

    // if my*Menu is currently posted, delet it/them now
    document.getElementById("myFileMenu"  ).classList.remove("showFileMenu"  );
    document.getElementById("myCapsMenu"  ).classList.remove("showCapsMenu"  );
    document.getElementById("myToolMenu"  ).classList.remove("showToolMenu"  );
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    if        (wv.curMode >= 8) {
        alert("stepThru not enabled when in tool");
    } else if (wv.curMode >= 0) {
        browserToServer("nextStep|"+direction+"|");

        var button = document.getElementById("stepThruBtn");
        button["innerHTML"] = "NextStep";
    } else {
        alert("Button disabled");
    }

    if (wv.curStep == 0) {
        postMessage("Entering StepThru: press \"n\" for next, \"p\" for previous, \"f\" for first, and \"l\" for last");
    }
};


//
// callback when "helpButton" is pressed (called by ESP.html)
//
var cmdHelp = function () {

    // if my*Menu is currently posted, delet it/them now
    document.getElementById("myFileMenu"  ).classList.remove("showFileMenu"  );
    document.getElementById("myCapsMenu"  ).classList.remove("showCapsMenu"  );
    document.getElementById("myToolMenu"  ).classList.remove("showToolMenu"  );
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    // open help in another tab
    window.open("ESP-help.html");
};


//
// callback when "testButton" is pressed (called by ESP.html)
//
var cmdTest = function () {
    // alert("in cmdTest()")

};


//
// callback when "collabButton" is pressed (called by ESP.html)
//
var cmdCollab = function () {
    // alert("in cmdCollab()");

    // toggle between hiding and showing the Collab menu contents
    document.getElementById("myCollabMenu").classList.toggle("showCollabMenu");

    // if myFileMenu or myCapsMenu is currently posted, delete it now
    document.getElementById("myFileMenu").classList.remove("showFileMenu");
    document.getElementById("myCapsMenu").classList.remove("showCapsMenu");

    // if myToolMenu is currently posted, delete it now
    document.getElementById("myToolMenu").classList.remove("showToolMenu");

    /* remove previous menu entries */
    var menu = document.getElementById("myCollabMenu");
    while (menu.firstChild !== null) {
        menu.removeChild(menu.firstChild)
    }

    /* add entries depending on collab mode */
    var button;

    // user with ball
    if (wv.myRole == 0) {
        var userList = wv.userNames.split("|");

        for (var i = 1; i < userList.length-1; i++) {
            if (userList[i] != wv.myName) {
                button = document.createElement("input");
                button.type    = "button";
                button.title   = "Pass";
                button.value   = "Pass ball to: "+userList[i];
                button.userindx= i;
                button.onclick = function (i) {cmdCollabPass(i);};
                menu.appendChild(button);
            }
        }

        button = document.createElement("input");
        button.type     = "button";
        button.title    = "Message";
        button.value    = "Send message";
        button.onclick  = cmdCollabMessage;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "List";
        button.value   = "List all users";
        button.onclick = cmdCollabListUsers;
        menu.appendChild(button);

    // sync'd user
    } else if (wv.myRole == 1) {
        button = document.createElement("input");
        button.type     = "button";
        button.title    = "Un-sync";
        button.value    = "Un-sync display";
        button.onclick  = cmdCollabUnsync;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type     = "button";
        button.title    = "Message";
        button.value    = "Send message";
        button.onclick  = cmdCollabMessage;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "List";
        button.value   = "List all users";
        button.onclick = cmdCollabListUsers;
        menu.appendChild(button);

    // un-sync'd user
    } else {
        button = document.createElement("input");
        button.type     = "button";
        button.title    = "Sync";
        button.value    = "Sync display";
        button.onclick  = cmdCollabSync;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type     = "button";
        button.title    = "Message";
        button.value    = "Send message";
        button.onclick  = cmdCollabMessage;
        menu.appendChild(button);

        button = document.createElement("input");
        button.type    = "button";
        button.title   = "List";
        button.value   = "List all users";
        button.onclick = cmdCollabListUsers;
        menu.appendChild(button);
    }

};


//
// callback when "Collab->List" is pressed
//
var cmdCollabListUsers = function () {

    // close the Collab menu
    var menu = document.getElementsByClassName("collabMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCollabMenu")) {
            menu[i].classList.remove(  "showCollabMenu");
        }
    }

    var mesg = "All users in this session:";

    var userList = wv.userNames.split("|");
    for (var i = 1; i < userList.length-1; i++) {
        mesg += "\n   "+userList[i];
    }

    if (wv.myRole == 0) {
        mesg += "\nI ("+wv.myName+") have the ball";
    } else if (wv.myRole == 1) {
        mesg += "\nI ("+wv.myName+") am sync'd";
    } else {
        mesg += "\nI ("+wv.myName+") am not sync'd";
    }

    alert(mesg);
};


//
// callback when "Collab->Pass" is pressed
//
var cmdCollabPass = function (e, indx) {
    // alert("in cmdCollabPass(e="+e+", indx="+indx+")");

    var index;
    if (e !== null) {
        index = e.srcElement.userindx;
    } else {
        index = indx;
    }

    // close the Collab menu
    var menu = document.getElementsByClassName("collabMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCollabMenu")) {
            menu[i].classList.remove(  "showCollabMenu");
        }
    }

    var userList = wv.userNames.split("|");

    // send message to server to pass the ball
    browserToServer("userName|"+wv.myName+"|"+userList[index]+"|");
};


//
// callback when "Collab->Sync" is pressed
//
var cmdCollabSync = function () {
    // alert("in cmdCollabSync()");

    // close the Collab menu
    var menu = document.getElementsByClassName("collabMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCollabMenu")) {
            menu[i].classList.remove(  "showCollabMenu");
        }
    }

    // change role
    wv.myRole = 1;

    // get last sent xform
    if (wv.lastXform !== null) {
        var readViewList = wv.lastXform.split("|");
        var entries      = readViewList[2].split(",");
        if (entries.length == 16) {
            var matrix  = Array(16);
            for (var i = 0; i < 16; i++) {
                matrix[i] = parseFloat(entries[i]);
            }
            wv.mvMatrix.makeIdentity();
            wv.uiMatrix.load(matrix);
            wv.scale    = parseFloat(readViewList[1]);
            wv.sceneUpd = 1;
        }
    }

    // ask user with ball to send its current state
    browserToServer("sendState|");

    // update button
    var button = document.getElementById("collabMenuBtn");
    button.style.backgroundColor = "#FFFFAF";    // yellow
};


//
// callback when "Collab->Unsync" is pressed
//
var cmdCollabUnsync = function () {
    // alert("in cmdCollabUnsync()");

    // close the Collab menu
    var menu = document.getElementsByClassName("collabMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCollabMenu")) {
            menu[i].classList.remove(  "showCollabMenu");
        }
    }

    wv.myRole = 2;

    var button = document.getElementById("collabMenuBtn");
    button.style.backgroundColor = null;       // white
};


//
// callback when "Collab->Message" is pressed
//
var cmdCollabMessage = function () {
    // alert("in cmdCollabMessage()");

    // close the Collab menu
    var menu = document.getElementsByClassName("collabMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showCollabMenu")) {
            menu[i].classList.remove(  "showCollabMenu");
        }
    }

    var message = prompt("Enter message to be sent to al users:");
    if (message !== null && message.length > 0) {
        browserToServer("message|"+wv.myName+" says: "+message+"|");
    }
};


//
// callback when "overlayBeg|timName|olayName|" message is received
//
var cmdOverlayBeg = function (text) {
    //alert("in cmdOverlayBeg(text="+text+")");

    var textList = text.split("|");

    if (wv.overlay !== undefined) {
        if        (textList[1] == "viewer"  && wv.overlay == viewer ) {
        } else if (textList[1] == "plotter" && wv.overlay == plotter) {
        } else {
            alert("cannot have recursive overlays");
        }
    } else {
        wv.timName = textList[0];

        if        (textList[1] == "viewer") {
            wv.overlay = viewer;
        } else if (textList[1] == "plotter") {
            wv.overlay = plotter;
        }

        // hide done menu button
        document.getElementById("doneMenuBtn").hidden = true;

        // hide the solver button
        document.getElementById("solveButton").hidden = true;

        // show the button that will be used to exit the overlay
        var button = document.getElementById("exitOlayBtn");
        button.hidden = false;
        button["innerHTML"] = "Exit " + textList[1];
        button.style.backgroundColor = "#3FFF3F";
    }
};


//
// callback when "ExitOverlay" is pressed
//
var cmdOverlayEnd = function () {
    // alert("in cmdOverlayEnd()");

    if (wv.overlay === undefined) {
        alert("there is no active overlay");
    } else if (checkIfWithBall() === false) {

    } else {

        // send the unlock message
        browserToServer("overlayEnd|"+wv.timName+"|");

        // unhide done menu button and hide overlay button
        document.getElementById("doneMenuBtn").hidden = false;
        document.getElementById("exitOlayBtn").hidden = true;

        // unhide the solver button
        document.getElementById("solveButton").hidden = false;

//$$$        alert("in cmdOverlayEnd: wv.overlay="+wv.overlay+"   wv.curTool.name="+wv.curTool.name);

        // we no longer have an overlay
        wv.overlay = undefined;
        wv.timName = "";
    }
};

////////////////////////////////////////////////////////////////////////////////


//
// callback when "Caps Value" is pressed in Tree
//
var addCval = function () {
    // alert("in addCval()");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the new name
    var name = prompt("Enter new Caps Value name:");
    if (name === null) {
        return;
    } else if (name.length <= 0) {
        return;
    }

    // check that name is valid
    if (name.match(/^[a-zA-Z][\w:]*$/) === null) {
        alert("'"+name+"' is not a valid name");
        return;
    }

    // check that the name does not exist already
    for (var icval = 0; icval < cval.length; icval++) {
        if (name == cval[icval].name) {
            alert("'"+name+"' already exists");
            return;
        }
    }

    // store the values locally
    var newCval = cval.length;

    cval[newCval] = new Array();

    cval[newCval].name = name;
    cval[newCval].nrow = 1;
    cval[newCval].ncol = 1;
    cval[newCval].value = new Array(1);

    cval[newCval].value[0] = "";

    // remember info for Parameter
    wv.curCval = newCval;

    // set up editValuForm
    if (setupEditValuForm() > 0) {
        return;
    }

    // hide the buttons associated with sensitivities
    document.getElementById("computeGeomSens").hidden = true;
    document.getElementById("computeTessSens").hidden = true;
    document.getElementById("clearVelocities").hidden = true;
    document.getElementById("deleteValue"    ).hidden = true;

    // post the editCval form (with the addValu header)
    changeMode(4);
};


//
// callback when Caps Value name is pressed in Tree
//
var editCval = function (e) {
    // alert("in editCval(e="+e+")");

    if        (wv.curMode == 5) {
        // currently editting another Parameter or Caps Value, so cancel (throwing away changes)
        editValuCancel();
    } else if (wv.curMode == 3) {
        // currently editting a Branch,          so cancel (throwing away changes)
        editBrchCancel();
    } else if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    wv.menuEvent = e;

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Caps Value name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    var jnode   = inode;
    var newName = name;
    while (myTree.parent[jnode] > 2) {
        jnode   = myTree.parent[jnode];
        newName = myTree.name[jnode] + newName;
    }
    name = newName.replace(/\u00a0/g, "");

    // get the Caps Value index
    var icval = -1;      // 0-bias
    var jcval;           // 1-bias (and temp)
    for (jcval = 0; jcval < cval.length; jcval++) {
        if (cval[jcval].name.replace(/\^/g, "") == name) {
            icval = jcval;
            break;
        }
    }

    if (icval < 0) {
        alert("|"+name+"| not found");
        return;
    } else {
        jcval = icval + 1;
    }

    // highlight this Parameter in the Tree
    var myElem = document.getElementById("node"+inode+"col1");
    myElem.className = "currentTD";

    // highlight all Branches that explicitly mention this Caps Value
    var re = RegExp("[\\][)(+\\-*/^.,;]"+cval[icval].name+"[\\][)(+\\-*/^.,;]");
    for (var ibrch = 0; ibrch < brch.length; ibrch++) {
        for (var iarg = 0; iarg < brch[ibrch].args.length; iarg++) {
            if (("("+brch[ibrch].args[iarg]+")").match(re) !== null) {
                for (jnode = 0; jnode < myTree.name.length; jnode++) {
                    var parent = myTree.parent[jnode];
                    while (parent > 0) {
                        if (parent == 3 &&
                            myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[ibrch].name) {
                            document.getElementById("node"+jnode+"col1").className = "childTD";
                        }
                        parent = myTree.parent[parent];
                    }
                }
            }
        }
    }

    // remember info for the current Parameter
    wv.curCval = icval;

    // set up editValuForm
    setupEditValuForm();

    // hide the buttons associated with sensitivities
    document.getElementById("computeGeomSens").hidden = true;
    document.getElementById("computeTessSens").hidden = true;
    document.getElementById("clearVelocities").hidden = true;
    document.getElementById("deleteValue"    ).hidden = true;

    // post the editCval form (with the editCval header)
    changeMode(5);
};


//
// callback when "Design Parameters" is pressed in Tree
//
var addPmtr = function () {
    // alert("in addPmtr()");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the new name
    var name = prompt("Enter new Parameter name");
    if (name === null) {
        return;
    } else if (name.length <= 0) {
        return;
    }

    // check that name is valid
    if (name.match(/^[a-zA-Z][\w:]*$/) === null) {
        alert("'"+name+"' is not a valid name");
        return;
    }

    // check that the name does not exist already
    for (var ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (name == pmtr[ipmtr].name) {
            alert("'"+name+"' already exists");
            return;
        }
    }

    // store the values locally
    var newPmtr = pmtr.length;

    pmtr[newPmtr] = new Array();

    pmtr[newPmtr].name = name;
    pmtr[newPmtr].type = OCSM_DESPMTR;
    pmtr[newPmtr].nrow = 1;
    pmtr[newPmtr].ncol = 1;
    pmtr[newPmtr].value = new Array(1);
    pmtr[newPmtr].dot   = new Array(1);

    pmtr[newPmtr].value[0] = "";
    pmtr[newPmtr].dot[  0] = "0";

    // remember info for Parameter
    wv.curPmtr = newPmtr;

    // set up editValuForm
    if (setupEditValuForm() > 0) {
        return;
    }

    // unhide the buttons associated with sensitivities
    document.getElementById("computeGeomSens").hidden = false;
    document.getElementById("computeTessSens").hidden = false;
    document.getElementById("clearVelocities").hidden = false;
    document.getElementById("deleteValue"    ).hidden = false;

    // post the editPmtr form (with the addPmtr header)
    changeMode(4);
};


//
// callback when Design Parameter name is pressed in Tree
//
var editPmtr = function (e) {
    // alert("in editPmtr(e="+e+")");

    if        (wv.curMode == 5) {
        // currently editting another Parameter, so cancel (throwing away changes)
        editValuCancel();
    } else if (wv.curMode == 3) {
        // currently editting a Branch,          so cancel (throwing away changes)
        editBrchCancel();
    } else if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    wv.menuEvent = e;

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Parameter name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    var jnode   = inode;
    var newName = name;
    while (myTree.parent[jnode] > 2) {
        jnode   = myTree.parent[jnode];
        newName = myTree.name[jnode] + newName;
    }
    name = newName.replace(/\u00a0/g, "");

    // get the Parameter index
    var ipmtr = -1;      // 0-bias
    var jpmtr;           // 1-bias (and temp)
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        if (pmtr[jpmtr].name.replace(/\^/g, "") == name) {
            ipmtr = jpmtr;
            break;
        }
    }

    if (ipmtr < 0) {
        alert("|"+name+"| not found");
        return;
    } else {
        jpmtr = ipmtr + 1;
    }

    // highlight this Parameter in the Tree
    var myElem = document.getElementById("node"+inode+"col1");
    myElem.className = "currentTD";

    // highlight all Branches that explicitly mention this Parameter
    var re = RegExp("[\\][)(+\\-*/^.,;]"+pmtr[ipmtr].name+"[\\][)(+\\-*/^.,;]");
    for (var ibrch = 0; ibrch < brch.length; ibrch++) {
        for (var iarg = 0; iarg < brch[ibrch].args.length; iarg++) {
            if (("("+brch[ibrch].args[iarg]+")").match(re) !== null) {
                for (jnode = 0; jnode < myTree.name.length; jnode++) {
                    var parent = myTree.parent[jnode];
                    while (parent > 0) {
                        if (parent == 3 &&
                            myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[ibrch].name) {
                            document.getElementById("node"+jnode+"col1").className = "childTD";
                        }
                        parent = myTree.parent[parent];
                    }
                }
            }
        }
    }

    // remember info for the current Parameter
    wv.curPmtr = ipmtr;

    // set up editValuForm
    setupEditValuForm();

    // unhide the buttons associated with sensitivities
    document.getElementById("computeGeomSens").hidden = false;
    document.getElementById("computeTessSens").hidden = false;
    document.getElementById("clearVelocities").hidden = false;
    document.getElementById("deleteValue"    ).hidden = false;

    // post the editPmtr form (with the editPmtr header)
    changeMode(5);
};


//
// callback when "Add row" is pressed in editValuForm (called by ESP.html)
//
var addRow = function () {
    // alert("in addRow()");

    // adding Paramater
    if (wv.curPmtr >= 0) {

        // adjust the number of rows
        pmtr[wv.curPmtr].nrow++;
        pmtr[wv.curPmtr].value = new Array(pmtr[wv.curPmtr].nrow*pmtr[wv.curPmtr].ncol);
        pmtr[wv.curPmtr].dot   = new Array(pmtr[wv.curPmtr].nrow*pmtr[wv.curPmtr].ncol);

        for (var i = 0; i < pmtr[wv.curPmtr].value.length; i++) {
            pmtr[wv.curPmtr].value[i] = "";
            pmtr[wv.curPmtr].dot[  i] = "0";
        }

    // adding Caps Value
    } else if (wv.curCval >= 0) {

        // adjust the number of rows
        cval[wv.curCval].nrow++;
        cval[wv.curCval].value = new Array(cval[wv.curCval].nrow*cval[wv.curCval].ncol);

        for (var i = 0; i < cval[wv.curCval].value.length; i++) {
            cval[wv.curCval].value[i] = "";
        }

    // ooops, problem
    } else {
        alert("we have a problem: curPmtr="+wv.curPmtr+" curCval="+wv.curCval);
    }

    // set up editValuForm
    if (setupEditValuForm() > 0) {
        return;
    }

    // post the editPmtr form (with the addPmtr header)
    changeMode(4);
};


//
// callback when "Add column" is pressed in editValuForm (called by ESP.html)
//
var addColumn = function () {
    // alert("in addColumn()");

    // adding Parameter
    if (wv.curPmtr >= 0) {

        // adjust the number of columns
        pmtr[wv.curPmtr].ncol++;
        pmtr[wv.curPmtr].value = new Array(pmtr[wv.curPmtr].nrow*pmtr[wv.curPmtr].ncol);
        pmtr[wv.curPmtr].dot   = new Array(pmtr[wv.curPmtr].nrow*pmtr[wv.curPmtr].ncol);

        for (var i = 0; i < pmtr[wv.curPmtr].value.length; i++) {
            pmtr[wv.curPmtr].value[i] = "";
            pmtr[wv.curPmtr].dot[  i] = "0";
        }

    // adding Caps Value
    } else if (wv.curCval >= 0) {

        // adjust the number of columns
        cval[wv.curCval].ncol++;
        cval[wv.curCval].value = new Array(cval[wv.curCval].nrow*cval[wv.curCval].ncol);

        for (var i = 0; i < cval[wv.curCval].value.length; i++) {
            cval[wv.curCval].value[i] = "";
        }

    // ooops, problem
    } else {
        alert("we have a problem: curPmtr="+wv.curPmtr+" curCval="+wv.curCval);
    }

    // set up editValuForm
    if (setupEditValuForm() > 0) {
        return;
    }

    // post the editPmtr form (with the addPmtr header)
    changeMode(4);
};


//
// callback when "Compute Geom Sens" is pressed in editValuForm (called by ESP.html)
//
var compGeomSens = function () {
    // alert("in compGeomSens()");

    if (checkIfWithBall() === false) {
        return;
    }

    // sensitivities cannot be computed for Caps Values
    if (wv.curCval >= 0) {
        alert("Sensitivities cannot be computed for a Caps Value");
        return;
    }

//    // disable this command if there were any changes to the Parameter
//    if (numberOfPmtrChanges() > 0) {
//        alert("Changes were made.  Press 'Cancel' or 'OK' first");
//        return;
//    }

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Parameter name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    var jnode   = inode;
    var newName = name;
    while (myTree.parent[jnode] > 4) {
        jnode   = myTree.parent[jnode];
        newName = myTree.name[jnode] + newName;
    }
    name = newName.replace(/\u00a0/g, "");

    // get the Parameter index
    var ipmtr = -1;      // 0-bias
    var jpmtr;           // 1-bias (and temp)
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        if (pmtr[jpmtr].name.replace(/\^/g, "") == name) {
            ipmtr = jpmtr;
            break;
        }
    }
    if (ipmtr < 0) {
        alert(name+" not found");
        return;
    } else {
        jpmtr = ipmtr + 1;
    }

    // can only compute sensitivity for a despmtr
    if (pmtr[ipmtr].type != OCSM_DESPMTR) {
        alert("Can only compute sensitivity for a DESPMTR (not a CFGPMTR).");
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // if a scalar...
    if (pmtr[ipmtr].nrow == 1 && pmtr[ipmtr].ncol == 1) {

        // clear any previous velocities
        browserToServer("clrVels|geom|");
        for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
            pmtr[jpmtr].dot[0] = 0;
        }

        // set velocity for ipmtr
        pmtr[ipmtr].dot[0] = 1.0;
        browserToServer("setVel|"+pmtr[ipmtr].name+"|1|1|1|");
        postMessage("Computing sensitivity with respect to "+pmtr[ipmtr].name);

    // if multi-valued, only compute if there is a non-zero velocity...
    } else {
        var editValuForm = document.getElementById("editValuForm");

        // count non-zero velocities
        var nonzero = 0;
        var irow;
        var icol;
        var index = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                var value   = myInput.value.replace(/\s/g, "");

                if (value.length <= 0) {
                    alert("Value at (row "+irow+", col "+icol+") is blank");
                    return;
                } else if (isNaN(value)) {
                    alert("Illegal number format in value at (row "+irow+", col "+icol+")");
                    return;
                } else if (Math.abs(Number(value)) > 1e-6) {
                    nonzero++;
                }
            }
        }
        if (nonzero == 0) {
            alert("Cannot compute sensitivity since no velocities (below) are non-zero");
            return;
        }

        // clear any previous velocities
        browserToServer("clrVels|geom|");
        for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
            pmtr[jpmtr].dot[0] = 0;
        }

        // set the non-zero velocities
        index = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                var value   = myInput.value.replace(/\s/g, "");

                browserToServer("setVel|"+pmtr[ipmtr].name+"|"+irow+"|"+icol+"|"+value+"|");
            }
        }
    }

    // rebuild
    browserToServer("build|0|");

    if (wv.server == "serveCAPS") {
        browserToServer("getCvals|7|");
        wv.cvalStat = 6000;
    }

    browserToServer("getPmtrs|");
    wv.pmtrStat = 6000;

    if (wv.server != "serveCAPS") {
        browserToServer("getBrchs|4|");
        wv.brchStat = 6000;
    }

    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode( 0);
    changeMode(-1);
};


//
// callback when "Compute Tess Sens" is pressed in editValuForm (called by ESP.html)
//
var compTessSens = function () {
    // alert("in compTessSens()");

    if (checkIfWithBall() === false) {
        return;
    }

    // sensitivities cannot be computed for Caps Values
    if (wv.curCval >= 0) {
        alert("Sensitivities cannot be computed for a Caps Value");
        return;
    }

//    // disable this command if there were any changes to the Parameter
//    if (numberOfPmtrChanges() > 0) {
//        alert("Changes were made.  Press 'Cancel' or 'OK' first");
//        return;
//    }

    // get the Tree Node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Parameter name
    var name = myTree.name[inode].replace(/\u00a0/g, "");
    name = name.replace(/\^/g, "");

    var jnode   = inode;
    var newName = name;
    while (myTree.parent[jnode] > 4) {
        jnode   = myTree.parent[jnode];
        newName = myTree.name[jnode] + newName;
    }
    name = newName.replace(/\u00a0/g, "");

    // get the Parameter index
    var ipmtr = -1;      // 0-bias
    var jpmtr;           // 1-bias (and temp)
    for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
        if (pmtr[jpmtr].name.replace(/\^/g, "") == name) {
            ipmtr = jpmtr;
            break;
        }
    }
    if (ipmtr < 0) {
        alert(name+" not found");
        return;
    } else {
        jpmtr = ipmtr + 1;
    }

    // can only compute sensitivity for a despmtr
    if (pmtr[ipmtr].type != OCSM_DESPMTR) {
        alert("Can only compute sensitivity for a DESPMTR (not a CFGPMTR).");
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // if a scalar...
    if (pmtr[ipmtr].nrow == 1 && pmtr[ipmtr].ncol == 1) {

        // clear any previous velocities
        browserToServer("clrVels|tess|");
        for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
            pmtr[jpmtr].dot[0] = 0;
        }

        // set velocity for ipmtr
        pmtr[ipmtr].dot[0] = 1.0;
        browserToServer("setVel|"+pmtr[ipmtr].name+"|1|1|1|");
        postMessage("Computing sensitivity with respect to "+pmtr[ipmtr].name);

    // if multi-valued, only compute if there is a non-zero velocity...
    } else {
        var editValuForm = document.getElementById("editValuForm");

        // count non-zero velocities
        var nonzero = 0;
        var irow;
        var icol;
        var index = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                var value   = myInput.value.replace(/\s/g, "");

                if (value.length <= 0) {
                    alert("Value at (row "+irow+", col "+icol+") is blank");
                    return;
                } else if (isNaN(value)) {
                    alert("Illegal number format in value at (row "+irow+", col "+icol+")");
                    return;
                } else if (Math.abs(Number(value)) > 1e-6) {
                    nonzero++;
                }
            }
        }
        if (nonzero == 0) {
            alert("Cannot compute sensitivity since no velocities (below) are non-zero");
            return;
        }

        // clear any previous velocities
        browserToServer("clrVels|tess|");
        for (jpmtr = 0; jpmtr < pmtr.length; jpmtr++) {
            pmtr[jpmtr].dot[0] = 0;
        }

        // set the non-zero velocities
        index = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                var value   = myInput.value.replace(/\s/g, "");

                browserToServer("setVel|"+pmtr[ipmtr].name+"|"+irow+"|"+icol+"|"+value+"|");
            }
        }
    }

    // rebuild
    browserToServer("build|0|");

    if (wv.server == "serveCAPS") {
        browserToServer("getCvals|");
        wv.cvalStat = 6000;
    }

    browserToServer("getPmtrs|");
    wv.pmtrStat = 6000;

    if (wv.server != "serveCAPS") {
        browserToServer("getBrchs|5|");
        wv.brchStat = 6000;
    }

    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode( 0);
    changeMode(-1);
};


//
// callback when "Set Design Velocity" is pressed in editValuForm (called by ESP.html)
//
var setVel = function () {
    // alert("in setVel()");

    alert("Set velocity(s) in table at bottom");

};


//
// callback when "Clear Design Velocities" is pressed in editValuForm (called by ESP.html)
//
var clrVels = function () {
    // alert("in clrVels()");

    if (checkIfWithBall() === false) {
        return;
    }

    // disable this command if there were any changes to the Parameter
    if (numberOfPmtrChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get an updated Parameter list (so that added Parameter is listed)
    browserToServer("clrVels|.|");

    // set all visible velocities to 0 and update display
    for (var ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (pmtr[ipmtr].type == OCSM_DESPMTR) {
            var index = -1;
            for (var irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
                for (var icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                    index++;

                    var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                    myInput.value = 0;

                    pmtr[ipmtr].dot[index] = 0;
                }
            }
        }
    }

    // update the UI
    postMessage("Design Velocities have been cleared ====> Re-build is needed <====");
    activateBuildButton();
};


//
// callback when "OK" is pressed in editValuForm (called by ESP.html)
//
var editValuOk = function () {
    // alert("in editValuOk()");

    if (checkIfWithBall() === false) {
        return;
    }

    var editValuForm = document.getElementById("editValuForm");

    // saving Parameter
    if (wv.curPmtr >= 0) {
        var ipmtr = wv.curPmtr;
        var name  = pmtr[ipmtr].name;
        var nrow  = pmtr[ipmtr].nrow;
        var ncol  = pmtr[ipmtr].ncol;
        var irow;
        var icol;

        // make sure that all entries have valid values
        var nchange = 0;

        var index   = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                // get the new value
                var myInput = editValuForm["row"+irow+"col"+icol+"val"];
                var value   = myInput.value.replace(/\s/g, "");

                if (value.length <= 0) {
                    alert("Value at (row "+irow+", col "+icol+") is blank");
                    return;
                } else if (isNaN(value)) {
                    alert("Illegal number format in value at (row "+irow+", col "+icol+")");
                    return;
                }
            }
        }

        if (pmtr[ipmtr].type == OCSM_DESPMTR) {
            index   = -1;
            for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                    index++;

                    // get the new value
                    var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                    var value   = myInput.value.replace(/\s/g, "");

                    if (value.length <= 0) {
                        alert("Velocity at (row "+irow+", col "+icol+") is blank");
                        return;
                    } else if (isNaN(value)) {
                        alert("Illegal number format in velocity at (row "+irow+", col "+icol+")");
                        return;
                    }
                }
            }
        }

        // send the new Parameter to the server if in add Pmtr mode
        if (wv.curMode == 4) {
            var mesg = "newPmtr|"+name+"|"+nrow+"|"+ncol+"|";

            index = -1;
            for (irow = 1; irow <= nrow; irow++) {
                for (icol = 1; icol <= ncol; icol++) {
                    index++;
                    mesg = mesg+"|";
                }
            }

            browserToServer(mesg);
        }

        // get each of the values
        index = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                // get the new value
                var myInput = editValuForm["row"+irow+"col"+icol+"val"];
                var value = myInput.value.replace(/\s/g, "");

                if (value != pmtr[ipmtr].value[index]) {
                    postMessage("Parameter '"+pmtr[ipmtr].name+"["+irow+","+icol+
                                "]' has been changed to "+value+
                                " ====> Re-build is needed <====");
                    nchange++;

                    // store the value locally
                    pmtr[ipmtr].value[index] = Number(value);

                    // send the new value to the server
                    browserToServer("setPmtr|"+pmtr[ipmtr].name+"|"+irow+"|"+icol+"|"+value+"|");
                }
            }
        }

        // get each of the velocities
        if (pmtr[ipmtr].type == OCSM_DESPMTR) {
            index = -1;
            for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
                for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                    index++;

                    // get the new value
                    var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                    var value = myInput.value.replace(/\s/g, "");

                    if (value != pmtr[ipmtr].dot[index]) {
                        postMessage("Velocity of parameter '"+pmtr[ipmtr].name+"["+irow+","+icol+
                                    "]' has been changed to "+value+
                                    " ====> Re-build is needed <====");
                        nchange++;

                        // store the value locally
                        pmtr[ipmtr].dot[index] = Number(value);

                        // send the new value to the server
                        browserToServer("setVel|"+pmtr[ipmtr].name+"|"+irow+"|"+icol+"|"+value+"|");
                    }
                }
            }
        }

        // update the UI
        if (nchange > 0) {
            wv.nchanges += nchange;

            if (wv.curMode != 4) {
                var id     = wv.menuEvent["target"].id;
                var myElem = document.getElementById(id);
                myElem.className = "fakelinkoff";

                // get an updated Parameter list (so that added Pmtr is listed)
            } else {
                browserToServer("getPmtrs|");
            }

            activateBuildButton();
        }

        wv.curPmtr = -1;

    // saving Caps Value
    } else if (wv.curCval >= 0) {
        var icval = wv.curCval;
        var name  = cval[icval].name;
        var nrow  = cval[icval].nrow;
        var ncol  = cval[icval].ncol;
        var irow;
        var icol;

        // make sure that all entries have valid values
        var nchange = 0;

        var index   = -1;
        for (irow = 1; irow <= cval[icval].nrow; irow++) {
            for (icol = 1; icol <= cval[icval].ncol; icol++) {
                index++;

                // get the new value
                var myInput = editValuForm["row"+irow+"col"+icol+"val"];
                var value   = myInput.value.replace(/\s/g, "");

                if (value.length <= 0) {
                    alert("Value at (row "+irow+", col "+icol+") is blank");
                    return;
                } else if (isNaN(value)) {
                    alert("Illegal number format in value at (row "+irow+", col "+icol+")");
                    return;
                }
            }
        }

        // send the new Parameter to the server if in add Cval mode
        if (wv.curMode == 4) {
            var mesg = "newCval|"+name+"|"+nrow+"|"+ncol+"|.|";

            index = -1;
            for (irow = 1; irow <= nrow; irow++) {
                for (icol = 1; icol <= ncol; icol++) {
                    index++;
                    mesg = mesg+"|";
                }
            }

            browserToServer(mesg);
        }

        // get each of the values
        index = -1;
        for (irow = 1; irow <= cval[icval].nrow; irow++) {
            for (icol = 1; icol <= cval[icval].ncol; icol++) {
                index++;

                // get the new value
                var myInput = editValuForm["row"+irow+"col"+icol+"val"];
                var value = myInput.value.replace(/\s/g, "");

                if (value != cval[icval].value[index]) {
                    nchange++;

                    // send the new value to the server
                    browserToServer("setCval|"+cval[icval].name+"|"+irow+"|"+icol+"|"+value+"|");
                }
            }
        }

        // update the UI
        if (nchange > 0) {
            wv.nchanges += nchange;

            if (wv.curMode != 4) {
                var id     = wv.menuEvent["target"].id;
                var myElem = document.getElementById(id);
                myElem.className = "fakelinkoff";

                // get an updated Caps Value list (so that added Cval is listed)
            } else {
                browserToServer("getCvals|");
            }
        }

        wv.curCval = -1;

    // ooops, problem
    } else {
        alert("we have a problem: curPmtr="+wv.curPmtr+" curCval="+wv.curCval);
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "Cancel" is pressed in editValuForm (called by ESP.html)
//
var editValuCancel = function () {
    // alert("in editValuCancel()");

    // if we are in process of adding a Parameter or Caps Value, remove it now
    if (wv.curMode == 4) {

        // canceling Parameter
        if (wv.curPmtr >= 0) {
            pmtr.splice(pmtr.length-1, 1);

            wv.curPmtr = -1;

        //canceling Caps Value
        } else if (wv.curCval >= 0) {
            cval.splice(cval.length-1, 1);

            wv.curCval = -1;

        // ooops, problem
        } else {
            alert("we have a problem: curPmtr="+wv.curPmtr+" curCval="+wv.curCval);
        }
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "Delete Parameter" or "Delete Caps Value" is pressed in editValuForm (called by ESP.html)
//
var delValu = function () {
    // alert("in delValu()");

    if (checkIfWithBall() === false) {
        return;
    }

    // deleting Parameter
    if (wv.curPmtr >= 0) {
        var ipmtr = wv.curPmtr + 1;

        // send message to the server
        browserToServer("delPmtr|"+pmtr[wv.curPmtr].name+"|");

        // get updated Parameters
        browserToServer("getPmtrs|");
        wv.pmtrStat = 0;

        // update the UI
        postMessage("Deleting Parameter "+name+" ====> Re-build is needed <====");
        activateBuildButton();

    // deleting Caps Value
    } else if (wv.curCval >= 0) {
        var icval = wv.curCval + 1;

        // send message to the server
        browserToServer("delCval|"+cval[wv.curCval].name+"|");

        // get updated Parameters
        browserToServer("getCvals|");
        wv.cvalStat = 0;

        // update the UI
        postMessage("Deleting Caps Value "+name+" ====> Re-build is needed <====");
        activateBuildButton();

    // ooops, problem
    } else {
        alert("we have a problem: curPmtr="+wv.curPmtr+" curCval="+wv.curCval);
    }

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "Local Variables" is pressed in Tree
//
var showOutpmtrs = function () {
    // alert("in showOutpmtrs()");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // count the number of OUTPMTRs
    var count = 0;
    for (var ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (pmtr[ipmtr].type == OCSM_OUTPMTR) {
            count++;
        }
    }

    if (count == 0) {
        alert("There are no OUTPMTRs");
        return;
    }

    // build up the message to display
    var message = "";
    for (ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (pmtr[ipmtr].type == OCSM_OUTPMTR) {
            message += pmtr[ipmtr].name + ":\n    Values:\n";
            var irow;
            var icol;
            var index = -1;
            for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
                message += "      ";
                for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                    index++;
                    message += " " + pmtr[ipmtr].value[index].toFixed(6).padStart(14);
                }
                message += "\n";
            }
            message += "    Velocities:\n";

            index = -1;
            for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
                message += "      ";
                for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                    index++;
                    message += " " + pmtr[ipmtr].dot[index].toFixed(6).padStart(14);
                }
                message += "\n";
            }
        }
    }

    // add the message to the outpmtrsList
    var outpmtrsList = document.getElementById("outpmtrsList");

    var pre  = document.createElement("pre");
    var text = document.createTextNode(message);
    pre.appendChild(text);
    outpmtrsList.insertBefore(pre, outpmtrsList.lastChild);

    // post the showOutpmtrsForm
    changeMode(6);
};


//
// callback when "Ok" is pressed in showOutpmtrsForm (called by ESP.html)
//
var showOutpmtrsOk = function () {
    // alert("in showOutpmtrsOk()");

    if (wv.curMode != 6) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    // clear the outpmtrsList
    var outpmtrsList = document.getElementById("outpmtrsList");
    while (outpmtrsList.firstChild !== null) {
        outpmtrsList.removeChild(outpmtrsList.firstChild)
    }

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "Branch" is pressed in Tree
// callback when "Add new Branch after this Brch" is pressed in editBrchForm (called by ESP.html)
//
var addBrch = function () {
    // alert("in addBrch()");

    // this check allows one to select AddBranchAfterThisBranch
    if (wv.menuEvent !== undefined && wv.curMode == 3) {
        if (brch[wv.curBrch].type == "udprim") {
            var arg0 = brch[wv.curBrch].args[0];
            if (arg0.charAt(1) == "$" || arg0.charAt(1) == "/") {
                alert("Cannot add a Branch within a UDC\nUse File->Edit instead.");
                return;
            }
        }
    } else if (wv.curMode != 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // remember
    if (wv.curBrch >= 0) {
        wv.afterBrch = wv.curBrch;
        wv.curBrch   = -1;
    } else {
        wv.afterBrch = brch.length - 1;
        wv.curBrch   = -1;
    }

    // post the addBrch form
    changeMode(1);
};


//
// callback when "OK" is pressed in addBrchForm (called by ESP.html)
//
var addBrchOk = function () {
    // alert("in addBrchOk()");

    // get the elements on the form
    var elements = document.getElementById("addBrchForm").elements;

    var jelem       = -1;
    for (var ielem = 0; ielem < elements.length; ielem++) {
        if (elements[ielem].checked) {
            jelem = ielem;
            break;
        }
    }

    // if nothing is selected, return to the WebViewer
    if (jelem < 0) {
        alert("Select a Branch type or press 'Cancel'");
        return;
    }

    // if something was picked, initialize the new Branch
    var newBrch = brch.length;

    brch[newBrch] = new Array();

    brch[newBrch].name = "**name_automatically_assigned**";
    brch[newBrch].type = elements[jelem].value;
    brch[newBrch].actv = OCSM_ACTIVE;
    brch[newBrch].args = new Array();

    // remember info for Branch to add after
    wv.curBrch = newBrch;

    // set up editBrchForm (and abort if problem is detected)
    if (setupEditBrchForm() > 0) {
        return;
    }

    // increment the number of changes
    wv.nchanges++;

    // post the editBrch form (with the addBrch header)
    changeMode(2);
};


//
// callback when "Cancel" is pressed in addBrchForm (called by ESP.html)
//
var addBrchCancel = function () {
    // alert("in addBrchCancel()");

    // return to the WebViewer
    changeMode(0);
};


//
// callback when Branch name is pressed in Tree
//
var editBrch = function (e) {
    // alert("in editBrch(e="+e+")");

    if        (wv.curMode == 3) {
        // currently editting another Branch, so cancel (throwing away changes)
        editBrchCancel();
    } else if (wv.curMode == 5) {
        // currently editting a Parameter,    so cancel (throwing away changes)
        editValuCancel();
    } else if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    changeMode(0);

    wv.menuEvent = e;

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    var name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    var ibrch = -1;           // 0-bias
    var jbrch;                // 1-bias (and temp)
    for (jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch;
            break;
        }
    }
    if (ibrch < 0) {
        alert(name+" not found");
        return;
    } else {
        jbrch = ibrch + 1;
    }

    // highlight this Branch in the Tree
    var myElem = document.getElementById("node"+inode+"col1");
    myElem.className = "currentTD";

    // highlight the parents and child of this Branch in the Tree */
    var ileft = brch[ibrch].ileft;
    var irite = brch[ibrch].irite;
    var ichld = brch[ibrch].ichld;

    for (var jnode = 0; jnode < myTree.name.length; jnode++) {
        if (ileft > 0) {
            if (// myTree.parent[jnode] == 3 &&        // why did I do this?
                myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[ileft-1].name) {
                myElem = document.getElementById("node"+jnode+"col1");
                myElem.className = "parentTD";
            }
        }
        if (irite > 0) {
            if (// myTree.parent[jnode] == 3 &&        // why did I do this?
                myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[irite-1].name) {
                myElem = document.getElementById("node"+jnode+"col1");
                myElem.className = "parentTD";
            }
        }
        if (ichld > 0) {
            if (// myTree.parent[jnode] == 3 &&        // why did I do this?
                myTree.name[jnode].replace(/\u00a0/g, "").replace(/>/g, "") == brch[ichld-1].name) {
                myElem = document.getElementById("node"+jnode+"col1");
                myElem.className = "childTD";
            }
        }
    }

    // remember info for the current Branch
    wv.curBrch   = ibrch;

    // set up editBrchForm
    if (setupEditBrchForm() > 0) {
        alert("This Branch cannot be editted\n(associated with a UDC)");
        unhighlightColumn1();
        changeMode(0);
        return;
    }

    // post the editBrch form (with the editBrch header)
    changeMode(3);
};


//
// callback when "Add Attribute" is pressed in editBrchForm  (called by ESP.html)
//
var addAttr = function () {
    // alert("in addAttr()");

    if (checkIfWithBall() === false) {
        return;
    }

    // disable this command if there were any chnges to the Branch
    if (numberOfBrchChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    var name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    var ibrch = -1;           // 0-bias
    var jbrch;                // 1-bias (and temp)
    for (jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch;
            break;
        }
    }
    if (ibrch < 0) {
        alert(name+" not found");
        return;
    } else {
        jbrch = ibrch + 1;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // get the Attribute name and value
    var atype  = prompt("Enter 1 for Attribute or 2 for Csystem", "1");
    if (atype === null) {
        return;
    } else if (atype == "1") {
        var aname  = prompt("Enter Attribute name");
        if (aname === null) {
            return;
        } else if (aname.length <= 0) {
            return;
        }
    } else if (atype == "2") {
        var aname  = prompt("Enter Csystem name");
        if (aname === null) {
            return;
        } else if (aname.length <= 0) {
            return;
        }
    } else {
        return;
    }
    var avalue = prompt("Enter Attribute value");
    if (avalue === null) {
        return;
    } else if (avalue.length <= 0) {
        return;
    }

    // send the new value to the server
    var mesg = "setAttr|"+jbrch+"|"+aname+"|"+atype+"|"+avalue+"|";

    browserToServer(mesg);

    // update the UI
    postMessage("Adding attribute '"+aname+"' (with value "+avalue+") to "+
                name+" ====> Re-build is needed <====");
    activateBuildButton();

    // get an updated version of the Branches
    wv.brchStat = 0;

    // return to  the WebViewer
    changeMode(0);
};


//
// callback when "Delete this Branch" is pressed in editBrchForm (called by ESP.html)
//
var delBrch = function () {
    // alert("in delBrch()");

    if (checkIfWithBall() === false) {
        return;
    }

    // disable this command if there were any changes to the Branch
    if (numberOfBrchChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    var ibrch = -1;
    var name  = "";

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // if this was called by pressing "Delete last Branch", adjust inode
    //    to point to last visible node (as if "Delete this Branch"
    //    was pressed on the last Branch)
    if (inode == 2) {
        inode = myTree.child[inode];

        while (myTree.next[inode] >= 0) {
            inode = myTree.next[inode];
        }
    }

    // get the Branch name
    name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    for (var jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch + 1;
        }
    }
    if (ibrch <= 0) {
        alert(name+" not found");
        return;
    }

    // do not allow a UDPARG or UDPRIM that points to a UDC to be deleted
    if (brch[ibrch-1].type == "udparg" || brch[ibrch-1].type == "udprim") {
        if (brch[ibrch-1].args[0].charAt(1) == "/" ||
            brch[ibrch-1].args[0].charAt(1) == "$"   ) {
            alert("Cannot delete \""+brch[ibrch-1].name+"\" that points to a UDC");
            return;
        }
    }

    // send the message to the server
    browserToServer("delBrch|"+ibrch+"|");

    // hide the Tree node that was just updated
    var element = myTree.document.getElementById("node"+inode);
    element.style.display = "none";

    // get an updated version of the Branches
    wv.brchStat = 0;

    // update the UI
    postMessage("Deleting Branch "+name+" ====> Re-build is needed <====");
    activateBuildButton();

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "Delete an Attribute" is pressed in editBrchForm (called by ESP.html)
//
var delAttr = function () {
    //alert("in delAttr");

    alert("To delete an Attribute, make its value blank");
};


//
// callback when "Show Attributes" is pressed in editBrchForm (called by ESP.html)
//
var showBrchAttrs = function () {
    // alert("in showBrchAttrs()");

    document.getElementById("AddBrchOrAttr").value   = "Add Attribute/Csystem";
    document.getElementById("AddBrchOrAttr").onclick = addAttr;

    document.getElementById("DelBrchOrAttr").value   = "Delete an Attribute/Csystem";
    document.getElementById("DelBrchOrAttr").onclick = delAttr;

    document.getElementById("ShowArgOrAttr").value   = "Show Arguments";
    document.getElementById("ShowArgOrAttr").onclick = showBrchArgs;

    document.getElementById("editBrchHeader2").innerHTML = "<h3>... or edit the attributes/csystems of the current Branch</h3>";

    document.getElementById("editBrchArgs" ).hidden = true;
    document.getElementById("editBrchAttrs").hidden = false;
};


//
// callback when "Show Arguments" is pressed in editBrchForm (called by ESP.html)
//
var showBrchArgs = function () {
    // alert("in showBrchArgs()");

    document.getElementById("AddBrchOrAttr").value   = "Add new Branch after this Branch";
    document.getElementById("AddBrchOrAttr").onclick = addBrch;

    document.getElementById("DelBrchOrAttr").value   = "Delete this Branch";
    document.getElementById("DelBrchOrAttr").onclick = delBrch;

    document.getElementById("ShowArgOrAttr").value   = "Show Attributes/Csystems";
    document.getElementById("ShowArgOrAttr").onclick = showBrchAttrs;

    document.getElementById("editBrchHeader2").innerHTML = "<h3>... or edit the arguments of the current Branch</h3>";

    document.getElementById("editBrchArgs" ).hidden = false;
    document.getElementById("editBrchAttrs").hidden = true;
};


//
// callback when "Build to the Branch" is pressed in editBrchForm (called by ESP.html)
//
var buildTo = function () {
    // alert("in buildTo()");

    if (checkIfWithBall() === false) {
        return;
    }

    // disable this command if there were any chnges to the Branch
    if (numberOfBrchChanges() > 0) {
        alert("Changes were made.  Press 'Cancel' or 'OK' first");
        return;
    }

    var ibrch = -1;
    var name  = "";

    // get the Tree node
    var id    = wv.menuEvent["target"].id;
    var inode = Number(id.substring(4,id.length-4));

    // get the Branch name
    name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

    // get the Branch index
    for (var jbrch = 0; jbrch < brch.length; jbrch++) {
        if (brch[jbrch].name == name) {
            ibrch = jbrch + 1;
        }
    }
    if (ibrch <= 0) {
        alert(name+" not found");
        return;
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // send the message to the server
    browserToServer("build|"+ibrch+"|");

    if (wv.server == "serveCAPS") {
        browserToServer("getCvals|8|");
        wv.cvalStat = 6000;
    }

    browserToServer("getPmtrs|");
    wv.pmtrStat = 6000;

    if (wv.server != "serveCAPS") {
        browserToServer("getBrchs|6|");
        wv.brchStat = 6000;
    }

    // update the UI
    postMessage("Re-building only to "+name+"...");

    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Re-building...";
    button.style.backgroundColor = "#FFFF3F";

    // inactivate buttons until build is done
    changeMode(-1);
};


//
// callback when "OK" button is pressed in editBrchForm (called by ESP.html)
//
var editBrchOk = function () {
    // alert("in editBrchOk()");

    if (checkIfWithBall() === false) {
        return;
    }

    var editBrchForm = document.getElementById("editBrchForm");

    var ibrch = wv.curBrch;
    var jbrch = ibrch + 1;

    var brchName = editBrchForm.brchName.value.replace(/\s/g, "");
    if (brchName.length <= 0) {
        alert("Name cannot be blank");
        return;
    }

    var mesg;
    var output;
    var nchange = 0;

    if (wv.curMode == 2) {
        mesg = "newBrch|"+(wv.afterBrch+1)+"|"+brch[ibrch].type+"|";
    } else {
        mesg = "setBrch|"+jbrch+"|";
    }

    // make sure that name does not contain a space or ">"
    var newBrchName = editBrchForm.brchName.value;

    if (newBrchName.indexOf(" ") >= 0) {
        alert("Changed name '"+newBrchName+"' cannot contain a space");
        return;
    } else if (newBrchName.indexOf(">") >= 0) {
        alert("Changed name '"+newBrchName+"' cannot contain '>'");
        return;
    }

//    // make sure that we are not adding or changing a Branch associated with a UDC
//    if (brch[ibrch].type == "udparg" || brch[ibrch].type == "udprim") {
//        var value = editBrchForm.argValu1.value.replace(/\s/g, "");
//        if (value.charAt(0) == "/" || value.charAt(0) == "$") {
//            if (wv.curMode == 2) {
//                alert("Cannot add a \""+brch[ibrch].type+"\" that calls a UDC\nUse File->Edit instead.");
//                return;
//            } else if (brch[ibrch].args[0] != "$"+value) {
//                alert("Cannot change primtype to a UDC\nUse File->Edit instead.");
//                return;
//            }
//        }
//    }

    if (newBrchName != brch[ibrch].name) {
        // make sure that name does not start with "Brch_"
        if (newBrchName.substring(0,5) == "Brch_") {
            alert("Changed name '"+newBrchName+
                  "' cannot begin with 'Brch_'");
            return;
        }

        // make sure that name does not already exist
        for (var kbrch = 0; kbrch < brch.length; kbrch++) {
            if (brch[kbrch].name == newBrchName) {
                alert("Name '"+newBrchName+"' already exists");
                return;
            }
        }

        // store the value locally
        postMessage("Changing name '"+brch[ibrch].name+"' to '"+
                    newBrchName+"' ====> Re-build is needed <====");
        brch[ibrch].name = newBrchName;

        nchange++;
    }

    if (wv.curMode == 3) {
        mesg = mesg + newBrchName + "|";
    }

    // update the activity
    var newActivity = editBrchForm.activity.value;
    if (brch[ibrch].actv != "none"     &&
        brch[ibrch].actv != newActivity  ) {

        // store the value locally
        brch[ibrch].actv = newActivity;
        if (brch[ibrch].actv == OCSM_ACTIVE) {
            postMessage("Activating Branch "+newBrchName+
                        " ====> Re-build is needed <====");
        } else {
            postMessage("Suppressing Branch "+newBrchName+
                        " ====> Re-build is needed <====");
        }

        // get an updated version of the Branches
        wv.brchStat = 0;
        nchange++;
    }

    if (wv.curMode == 2) {

    } else if (newActivity == OCSM_ACTIVE) {
        mesg = mesg + "active|";
    } else {
        mesg = mesg + "suppressed|";
    }

    // update any necessary arguments
    var prev_value = "dum";
    for (var iarg = 0; iarg < wv.numArgs; iarg++) {
        if        (iarg == 0) {
            var name  = document.getElementById("argName1").firstChild["data"];
            var value = editBrchForm.            argValu1.value.replace(/\s/g, "");
        } else if (iarg == 1) {
            var name  = document.getElementById("argName2").firstChild["data"];
            var value = editBrchForm.            argValu2.value.replace(/\s/g, "");
        } else if (iarg == 2) {
            var name  = document.getElementById("argName3").firstChild["data"];
            var value = editBrchForm.            argValu3.value.replace(/\s/g, "");
        } else if (iarg == 3) {
            var name  = document.getElementById("argName4").firstChild["data"];
            var value = editBrchForm.            argValu4.value.replace(/\s/g, "");
        } else if (iarg == 4) {
            var name  = document.getElementById("argName5").firstChild["data"];
            var value = editBrchForm.            argValu5.value.replace(/\s/g, "");
        } else if (iarg == 5) {
            var name  = document.getElementById("argName6").firstChild["data"];
            var value = editBrchForm.            argValu6.value.replace(/\s/g, "");
        } else if (iarg == 6) {
            var name  = document.getElementById("argName7").firstChild["data"];
            var value = editBrchForm.            argValu7.value.replace(/\s/g, "");
        } else if (iarg == 7) {
            var name  = document.getElementById("argName8").firstChild["data"];
            var value = editBrchForm.            argValu8.value.replace(/\s/g, "");
        } else if (iarg == 8) {
            var name  = document.getElementById("argName9").firstChild["data"];
            var value = editBrchForm.            argValu9.value.replace(/\s/g, "");
        }

        // make sure non-blank does not follow blank
        if (value.length > 0 && prev_value.length == 0) {
            alert(name+" should be blank (follows blank)");
            return;
        }

        // check for blanks (allowed for some Branch types because they
        // have variable number of arguments)
        if (value.length <= 0) {
            if (brch[ibrch].type == "evaluate") {
                if (iarg == 0) {
                    alert(name+" should not be blank");
                    return;
                }
            } else if (brch[ibrch].type == "udparg" || brch[ibrch].type == "udprim") {
                if (iarg == 0) {
                    alert(name+" should not be blank");
                    return;
                } else if (iarg == 2 || iarg == 4 || iarg == 6 || iarg == 8) {
                    if (prev_value.length > 0) {
                        alert(name+" should not be blank (follows non-blank)");
                        return;
                    }
                }
            } else if (brch[ibrch].type == "select") {
                if (iarg == 0) {
                    alert(name+" should not be blank");
                    return;
                }
            } else {
                alert(name+" should not be blank");
                return;
            }
        }

        // prepend a dollar sign if a string argument (ie, the name starts with a dollar sign)
        if (name.charAt(0) != "$" || value.length <= 0) {
            output =       value;
        } else {
            output = "$" + value;
        }

        // check for bars
        if (value.indexOf("|") >= 0) {
            alert(name+" should not contain a bar");
            return;
        }

        // check if there were any changes
        if (output != brch[ibrch].args[iarg]) {
            if (wv.curMode == 3) {
                postMessage("Changing "+name+" from '"+brch[ibrch].args[iarg]+
                            "' to '"+output+
                            "' ====> Re-build is needed <====");
            }

            brch[ibrch].args[iarg] = output;
            nchange++;
        }

        mesg = mesg + output + "|";

        // save previous value
        prev_value = value;
    }

    // put rest of bars at end of mesg
    for (iarg = wv.numArgs; iarg < 9; iarg++) {
        mesg = mesg + "|";
    }

    if (wv.curMode == 3) {
        mesg = mesg + "||";
    }

    // if adding a new Branch, send the "newBrch" message
    if (wv.curMode == 2) {
        postMessage("Branch (type="+brch[ibrch].type+") has been added"+
                "  ====> Re-build is needed <====");

        browserToServer(mesg);

    // if there are changes, sent the "setBrch" message
    } else if (nchange > 0) {
        browserToServer(mesg);
    }

    // update any necessary attributes
    if (wv.curMode != 2) {
        for (var iattr = 0; iattr < brch[ibrch].attrs.length; iattr++) {
            if        (iattr == 0) {
                var name  = document.getElementById("attrName1").firstChild["data"];
                var value = editBrchForm.            attrValu1.value.replace(/\s/g, "");
            } else if (iattr == 1) {
                var name  = document.getElementById("attrName2").firstChild["data"];
                var value = editBrchForm.            attrValu2.value.replace(/\s/g, "");
            } else if (iattr == 2) {
                var name  = document.getElementById("attrName3").firstChild["data"];
                var value = editBrchForm.            attrValu3.value.replace(/\s/g, "");
            } else if (iattr == 3) {
                var name  = document.getElementById("attrName4").firstChild["data"];
                var value = editBrchForm.            attrValu4.value.replace(/\s/g, "");
            } else if (iattr == 4) {
                var name  = document.getElementById("attrName5").firstChild["data"];
                var value = editBrchForm.            attrValu5.value.replace(/\s/g, "");
            } else if (iattr == 5) {
                var name  = document.getElementById("attrName6").firstChild["data"];
                var value = editBrchForm.            attrValu6.value.replace(/\s/g, "");
            } else if (iattr == 6) {
                var name  = document.getElementById("attrName7").firstChild["data"];
                var value = editBrchForm.            attrValu7.value.replace(/\s/g, "");
            } else if (iattr == 7) {
                var name  = document.getElementById("attrName8").firstChild["data"];
                var value = editBrchForm.            attrValu8.value.replace(/\s/g, "");
            } else if (iattr == 8) {
                var name  = document.getElementById("attrName9").firstChild["data"];
                var value = editBrchForm.            attrValu9.value.replace(/\s/g, "");
            }

            // check for bars
            if (value.indexOf("|") >= 0) {
                alert(name+" should not contain a bar");
                return;
            }

            // check if there were any changes
            if (value != brch[ibrch].attrs[iattr][2]) {
                if        (brch[ibrch].attrs[iattr][1] == "(attr)") {
                    if (value.length > 0) {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|1|"+value+"|");
                        postMessage("Changing attribute '"+name+"' to "+
                                    value+" ====> Re-build is needed <====");
                    } else {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|1|<DeLeTe>|");
                        postMessage("Attribute '"+name+
                                    "' being deleted. ====> Re-build is needed <====");
                    }
                } else if (brch[ibrch].attrs[iattr][1] == "(csys)") {
                    if (value.length > 0) {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|2|"+value+"|");
                        postMessage("Changing csystem '"+name+"' to "+
                                    value+" ====> Re-build is needed <====");
                    } else {
                        browserToServer("setAttr|"+jbrch+"|"+name+"|2|<DeLeTe>|");
                        postMessage("Csystem '"+name+
                                    "' being deleted. ====> Re-build is needed <====");
                    }
                } else {
                    alert("ERROR:: type="+brch[ibrch].attrs[iattr][1]);
                }

                brch[ibrch].attrs[iattr][2] = value;
                nchange++;
            }
        }
    }

    if (nchange > 0 || wv.curMode == 2) {
        wv.nchanges += nchange;

        // get an updated Branch list (so that added Branch is listed)
        browserToServer("getBrchs|7|");

        // update the UI
        if (wv.curMode == 3) {
            var id     = wv.menuEvent["target"].id;
            var myElem = document.getElementById(id);
            myElem.className = "fakelinkoff";
        }

        activateBuildButton();
    } else {
        alert("no changes were made");
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // if not adding skbeg/skend, return to the WebViewer
    if (wv.curMode != 2 || brch[ibrch].type != "skbeg") {
        changeMode(0);

    //otherwise enter the sketcher
    } else {
        sketch.ibrch = ibrch + 1;

        // send the message to the server
        browserToServer("loadSketch|"+sketch.ibrch+"|");

        // inactivate buttons until build is done
        changeMode(-1);
    }
};


//
// callback when "Cancel" is pressed in editBrchForm (called by ESP.html)
//
var editBrchCancel = function () {
    // alert("in editBrchCancel()");

    // if we are in process of adding a Branch, remove it now
    if (wv.curMode == 2) {
        brch.splice(brch.length-1, 1);
    }

    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // return to the WebViewer
    changeMode(0);
};


//
// callback when "ExpandAll" or "CollapseAll" is pressed for CapsValues in Tree
//
var allCvals = function () {
    // alert("in allCvals()");

    for (var jnode = 1; jnode < myTree.name.length; jnode++) {
        if (myTree.name[jnode] == "Caps Values") {
            if (myTree.prop1[jnode] == "ExpandAll") {
                for (var inode = 0; inode < myTree.parent.length-1; inode++) {
                    var parent = myTree.parent[inode];
                    while (parent > 0) {
                        if (parent == jnode) {
                            myTree.opened[inode] = 1;
                            break;
                        }
                        parent = myTree.parent[parent];
                    }
                }

                myTree.opened[jnode] = 1;
                myTree.prop1[ jnode] = "CollapseAll";

                var thisNode = document.getElementById("node"+jnode+"col3")
                thisNode["innerHTML"] = myTree.prop1[jnode];

                myTree.update();
            } else {
                //      myTree.prop1[jnode] = "ExpandAll";   // automagically done by rebuildTreeWindow();
                myTree.opened[jnode] = 0;
                rebuildTreeWindow(7);
            }
            break;
        }
    }
};


//
// callback when "ExpandAll" or "CollapseAll" is pressed for DesignParameters in Tree
//
var allPmtrs = function () {
    // alert("in allPmtrs()");

    for (var jnode = 1; jnode < myTree.name.length; jnode++) {
        if (myTree.name[jnode] == "Design Parameters") {
            if (myTree.prop1[jnode] == "ExpandAll") {
                for (var inode = 0; inode < myTree.parent.length-1; inode++) {
                    var parent = myTree.parent[inode];
                    while (parent > 0) {
                        if (parent == jnode) {
                            myTree.opened[inode] = 1;
                            break;
                        }
                        parent = myTree.parent[parent];
                    }
                }

                myTree.opened[jnode] = 1;
                myTree.prop1[ jnode] = "CollapseAll";

                var thisNode = document.getElementById("node"+jnode+"col3")
                thisNode["innerHTML"] = myTree.prop1[jnode];

                myTree.update();
            } else {
                //      myTree.prop1[jnode] = "ExpandAll";   // automagically done by rebuildTreeWindow();
                myTree.opened[jnode] = 0;
                rebuildTreeWindow(8);
            }
            break;
        }
    }
};


//
// callback when "ExpandAll" or "CollapseAll" is pressed for LocalVariables in Tree
//
var allLocals = function () {
    // alert("in allLocals()");

    for (var jnode = 1; jnode < myTree.name.length; jnode++) {
        if (myTree.name[jnode] == "Local Variables") {
            if (myTree.prop1[jnode] == "ExpandAll") {
                for (var inode = 0; inode < myTree.parent.length-1; inode++) {
                    var parent = myTree.parent[inode];
                    while (parent > 0) {
                        if (parent == jnode) {
                            myTree.opened[inode] = 1;
                            break;
                        }
                        parent = myTree.parent[parent];
                    }
                }

                myTree.opened[jnode] = 1;
                myTree.prop1[ jnode] = "CollapseAll";

                var thisNode = document.getElementById("node"+jnode+"col3")
                thisNode["innerHTML"] = myTree.prop1[jnode];

                myTree.update();
            } else {
                //      myTree.prop1[jnode] = "ExpandAll";   // automagically done by rebuildTreeWindow();
                myTree.opened[jnode] = 0;
                rebuildTreeWindow(9);
            }
            break;
        }
    }
};


//
// callback when "ExpandAll" or "CollapseAll" is pressed for Branches in Tree
//
var allBrchs = function () {
    // alert("in allBrchs()");

    for (var jnode = 1; jnode < myTree.name.length; jnode++) {
        if (myTree.name[jnode] == "Branches") {
            if (myTree.prop1[jnode] == "ExpandAll") {
                for (var inode = 0; inode < myTree.parent.length-1; inode++) {
                    var parent = myTree.parent[inode];
                    while (parent > 0) {
                        if (parent == jnode) {
                            myTree.opened[inode] = 1;
                            break;
                        }
                        parent = myTree.parent[parent];
                    }
                }

                myTree.opened[jnode] = 1;
                myTree.prop1[ jnode] = "CollapseAll";

                var thisNode = document.getElementById("node"+jnode+"col3")
                thisNode["innerHTML"] = myTree.prop1[jnode];

                myTree.update();
            } else {
                //      myTree.prop1[jnode] = "ExpandAll";   // automagically done by rebuildTreeWindow();
                myTree.opened[jnode] = 0;
                rebuildTreeWindow(10);
            }
            break;
        }
    }
};


//
// callback when "Display" is pressed in Tree
//
var chgDisplay = function () {
    // alert("in chgDislay()");

    if (checkIfFree() === false) {
        return;
    }

    var change = prompt("Enter type of display change:\n"+
                        "  +1 show Nodes\n"+
                        "  -1 hide Nodes\n"+
                        "  +2 show Edges\n"+
                        "  -2 hide Edges\n"+
                        "  +3 show Faces\n"+
                        "  -3 hide Faces\n"+
                        "  +4 show Csystems\n"+
                        "  -4 hide Csystems", "0");

    if (change == null) {
        return;
    } else if (isNaN(change)) {
        alert("Illegal number format, so no change");
        return;
    } else if (Number(change) == +1) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Nodes/)) {
                changeProp(inode, 1, "on");
            }
        }

        myTree.update();
    } else if (Number(change) == -1) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Nodes/)) {
                changeProp(inode, 1, "off");
            }
        }

        myTree.update();
    } else if (Number(change) == +2) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Edges/)) {
                changeProp(inode, 1, "on");
            }
        }

        myTree.update();
    } else if (Number(change) == -2) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Edges/)) {
                changeProp(inode, 1, "off");
            }
        }

        myTree.update();
    } else if (Number(change) == +3) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Faces/)) {
                changeProp(inode, 1, "on");
            }
        }

        myTree.update();
    } else if (Number(change) == -3) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Faces/)) {
                changeProp(inode, 1, "off");
            }
        }

        myTree.update();
    } else if (Number(change) == +4) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Csystems/)) {
                changeProp(inode, 1, "on");
            }
        }

        myTree.update();
    } else if (Number(change) == -4) {
        for (var inode = 0; inode < myTree.parent.length; inode++) {
            var name = myTree.name[myTree.parent[inode]];
            if (name !== undefined && name.match(/\u00a0*Csystems/)) {
                changeProp(inode, 1, "off");
            }
        }

        myTree.update();
    }
}


//
// callback when Body name is pressed in Tree
//
var showBodyAttrs = function (e) {
    // alert("in showBodyAttrs()");

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    var myElem = document.getElementById("node"+inode+"col2");
    var bodyName = myElem.firstChild.nodeValue;
    while (bodyName.charCodeAt(0) == 160) {
        bodyName = bodyName.slice(1);
    }

    var mesg  = bodyName+":";
    try {
        var attrs = wv.sgData[bodyName];
        for (var i = 0; i < attrs.length; i+=2) {
            mesg = mesg + "\n        "+attrs[i]+"= "+attrs[i+1];
        }
    } catch (x) {
    }
    postMessage(mesg);
};

////////////////////////////////////////////////////////////////////////////////

main.name = "main";


//
// callback when "homeButton" is pressed in main mode
//
main.cmdHome = function() {
    wv.mvMatrix.makeIdentity();
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when "leftButton" is pressed in main mode
//
main.cmdLeft = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(+90, 0,1,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when "riteButton" is pressed in main mode
//
main.cmdRite = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(-90, 0,1,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when "botmButton" is pressed in main mode
//
main.cmdBotm = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(-90, 1,0,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when "topButton" is pressed in main mode
//
main.cmdTop = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(+90, 1,0,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when "inButton" is pressed in main mode
//
main.cmdIn = function() {
    wv.mvMatrix.load(wv.uiMatrix);
    wv.mvMatrix.scale(2.0, 2.0, 2.0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale   *= 2.0;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when "outButton" is pressed in main mode
//
main.cmdOut = function() {
    wv.mvMatrix.load(wv.uiMatrix);
    wv.mvMatrix.scale(0.5, 0.5, 0.5);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale   *= 0.5;
    wv.sceneUpd = 1;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when any mouse is pressed in canvas in main mode
//
main.mouseDown = function(e) {
    if (!e) var e = event;

    if (wv.myRole == 1) return;

    wv.startX   =  e.clientX - wv.offLeft             - 1;
    wv.startY   = -e.clientY + wv.offTop  + wv.height + 1;

    wv.dragging = true;
    wv.button   = e.button;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
};


//
// callback when the mouse moves in canvas in main mode
//
main.mouseMove = function(e) {
    if (!e) var e = event;

    if (wv.myRole == 1) return;

    wv.cursorX  =  e.clientX - wv.offLeft             - 1;
    wv.cursorY  = -e.clientY + wv.offTop  + wv.height + 1;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;
};


//
// callback when the mouse is released in canvas in main mode
//
main.mouseUp = function(e) {
    wv.dragging = false;

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when the mouse wheel is rolled in canvas in main mode
//
main.mouseWheel = function(e) {
    if (e) {

        if (wv.myRole == 1) return;

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

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when the mouse leaves the canvas in main mode
//
main.mouseLeftCanvas = function(e) {
    if (wv.dragging) {
        wv.dragging = false;
    }

    // send matrix to server in case there are any users who are sync'd
    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("xform|"+wv.width+"|"+wv.height+"|"+wv.scale+"|"+wv.uiMatrix.getAsArray()+"|");
    }
};


//
// callback when a key is pressed in main mode
//
main.keyPress = function(e) {
//    if (!e) var e = event;

    // if <esc> was pressed, return to base mode
    if (e.charCode == 0 && e.keyCode == 27) {
        activateBuildButton();
        changeMode(0);
    }

    // if in canvas, record info about keypress
    if (wv.usingMain == 1) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

                        wv.modifier  = 0;
        if (e.shiftKey) wv.modifier |= 1;
        if (e.altKey  ) wv.modifier |= 2;
        if (e.ctrlKey ) wv.modifier |= 4;

    // if addBrchForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 1) {
        if (e.keyCode == 13) {
            addBrchOk();
            return false;
        }

    // if editBrchForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 2 || wv.curMode == 3) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

        if (wv.keyCode == 13) {
            editBrchOk();
            return false;
        } else {
            return true;
        }

    // if editValuForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 4 || wv.curMode == 5) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

        if (wv.keyCode == 13) {
            editValuOk();
            return false;
        } else {
            return true;
        }

    // if showOutpmtrsForm is posted, press OK when <return> is pressed
    } else if (wv.curMode == 6) {
        wv.keyPress = e.charCode;
        wv.keyCode  = e.keyCode;

        if (wv.keyCode == 13) {
            showOutpmtrsOk();
            return false;
        } else {
            return true;
        }
    }

    return true;
};


//
// callback when an arrow... or shift is pressed (needed for Chrome)
//
main.keyDown = function(e) {
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
};


//
// callback when a shift is released (needed for Chrome)
//
main.keyUp = function(e) {

    if (e.charCode == 0 && e.keyCode == 16) {
        wv.modifier = 0;
    }
};


//
// callback when "solveButton" is pressed
//
main.cmdSolve = function () {
    // alert("in cmdSolve()");

    var button  = document.getElementById("solveButton");
    var buttext = button["innerHTML"];

    if (buttext == "Up to date") {
        if (confirm("The configuration is up to date.\n" +
                    "Do you want to force a rebuild?") === true) {
            postMessage("Forced re-building...");

            // build first so that parameters are updated
            browserToServer("build|-1|");

            if (wv.server == "serveCAPS") {
                browserToServer("getCvals|9|");
                wv.cvalStat = 6000;
            }

            browserToServer("getPmtrs|");
            wv.pmtrStat = 6000;

            if (wv.server != "serveCAPS") {
                browserToServer("getBrchs|8|");
                wv.brchStat = 6000;
            }

            button["innerHTML"] = "Re-building...";
            button.style.backgroundColor = "#FFFF3F";

            // turn the background of the message window back to original color
            var botm = document.getElementById("brframe");
            botm.style.backgroundColor = "#F7F7F7";

            // inactivate buttons until build is done
            changeMode(-1);
        }

    } else if (buttext == "Press to Re-build") {
        if (wv.curMode != 0) {
            alert("Command disabled,  Press 'Cancel' or 'OK' first");
            return;
        }

        // build first so that parameters are updated
        browserToServer("build|0|");

        if (wv.server == "serveCAPS") {
            browserToServer("getCvals|10|");
            wv.cvalStat = 6000;
        }

        browserToServer("getPmtrs|");
        wv.pmtrStat = 6000;

        if (wv.server != "serveCAPS") {
            browserToServer("getBrchs|9|");
            wv.brchStat = 6000;
        }

        button["innerHTML"] = "Re-building...";
        button.style.backgroundColor = "#FFFF3F";

        // turn the background of the message window back to original color
        var botm = document.getElementById("brframe");
        botm.style.backgroundColor = "#F7F7F7";

        // inactivate buttons until build is done
        changeMode(-1);

    } else if (buttext == "Re-building...") {
        if (confirm("A rebuild is in progress.\n" +
                   "Do you really want to reset the ESP interface?") === true) {
            postMessage("WARNING:: User requested that ESP interface be reset...");

            activateBuildButton();
            changeMode(0);
        }

    } else if (buttext == "Fix before re-build") {
        alert("Edit/add/delete a Branch to fix before re-building");

    } else {
        alert("Unexpected button text:\""+buttext+"\"");
    }
};


//
// callback when "undoButton" is pressed (called by ESP.html)
//
main.cmdUndo = function () {
    // alert("in cmdUndo()");

    if (checkIfWithBall() === false) {
        return;
    }

    // if my*Menu is currently posted, delet it/them now
    document.getElementById("myFileMenu"  ).classList.remove("showFileMenu"  );
    document.getElementById("myToolMenu"  ).classList.remove("showToolMenu"  );
    document.getElementById("myCollabMenu").classList.remove("showCollabMenu");

    if (wv.curMode != 0) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    }

    browserToServer("undo|");
};


//
// callback when the mouse is pressed in key window
//
var setKeyLimits = function (e) {
    if (wv.curMode != 0) {
        alert("KeyLimits can only be currrently set in main mode");
        return;
    }

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

    if (Number(tempup) <= Number(templo)) {
        alert("Upper limit must be greater than lower limit");
        return;
    }

    if (templo != wv.loLimit || tempup != wv.uplimit) {
        wv.loLimit = templo;
        wv.upLimit = tempup;

        // send the limits back to the server
        browserToServer("setLims|"+wv.plotType+"|"+wv.loLimit+"|"+wv.upLimit+"|");
    }
};


//
// callback when the mouse is double-clicked in messages window
//
var gotoCsmError = function (e) {

    var botm = document.getElementById("brframe");

    // look for [[filename:linenum]] at the cursor
    var thisLine = e.rangeParent.data;
    var beg = thisLine.substring(0,e.rangeOffset).lastIndexOf("[[");
    var end = thisLine.indexOf("]]", e.rangeOffset);

    if (beg >= 0 && end > beg) {
        var foo = thisLine.slice(beg+2, end).split(":");
        if (foo.length == 2) {
            var filelist = wv.filenames.split("|");
            for (var ielem = 0; ielem < filelist.length; ielem++) {
                if (filelist[ielem] == foo[0]) {
                    wv.linenum = Number(foo[1]);
                    cmdFileEdit(null, ielem);
                    return;
                }
            }
        }
    }

    // if not found, look for last [[filename:lnenum]]
    var msgText = botm.innerText;
    beg = msgText.lastIndexOf("[[");
    end = msgText.lastIndexOf("]]");

    if (beg >= 0 && end > beg) {
        foo = msgText.slice(beg+2, end).split(":");
        if (foo.length == 2) {
            filelist = wv.filenames.split("|");
            for (var ielem = 0; ielem < filelist.length; ielem++) {
                if (filelist[ielem] == foo[0]) {
                    wv.linenum = Number(foo[1]);
                    cmdFileEdit(null, ielem);
                    return;
                }
            }
        } else {
            postMessage("\""+foo+"\" could not be parsed");
        }
    } else {
        postMessage("no [[filename:linenum]] found");
    }
};


//
// callback to toggle Viz property
//
var toggleViz = function (e) {
    // alert("in toggleViz(e="+e+")");

    if (wv.usingMain != 1) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (checkIfFree() === false) {
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Viz property
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
            changeProp(inode, 1, "on");
        } else {
            changeProp(inode, 1, "off");
        }

    //  toggle the Viz property (on all Faces/Edges/Nodes in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col3");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            changeProp(inode, 1, "on");
        } else {
            changeProp(inode, 1, "off");
        }
    }

    myTree.update();
};


//
// callback to toggle Grd property
//
var toggleGrd = function (e) {
    // alert("in toggleGrd(e="+e+")");

    if (wv.usingMain != 1) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (checkIfFree() === false) {
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Grd property
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
            changeProp(inode, 2, "on");
        } else {
            changeProp(inode, 2, "off");
        }

    // toggle the Grd property (on all Faces/Edges/Nodes in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col4");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            changeProp(inode, 2, "on");
        } else {
            changeProp(inode, 2, "off");
        }
    }

    myTree.update();
};


//
// callback to toggle Trn property
//
var toggleTrn = function (e) {
    // alert("in toggleTrn(e="+e+")");

    if (wv.usingMain != 1) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (checkIfFree() === false) {
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Trn property (on a Face)
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
            changeProp(inode, 3, "on");
        } else {
            changeProp(inode, 3, "off");
        }

    // toggle the Trn property (on all Faces in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col5");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            changeProp(inode, 3, "on");
            myElem.setAttribute("class", "fakelinkon");
            myElem.title = "Toggle Transparency off";
        } else {
            changeProp(inode, 3, "off");
            myElem.setAttribute("class", "fakelinkoff");
            myElem.title = "Toggle Transparency on";
        }
    }

    myTree.update();
};


//
// callback to toggle Ori property
//
var toggleOri = function (e) {
    // alert("in toggleOri(e="+e+")");

    if (wv.usingMain != 1) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (checkIfFree() === false) {
        return;
    }

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Ori property (on an Edge)
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
            changeProp(inode, 3, "on");
        } else {
            changeProp(inode, 3, "off");
        }

    // toggle the Ori property (on all Edges in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col5");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            changeProp(inode, 3, "on");
            myElem.setAttribute("class", "fakelinkon");
            myElem.title = "Toggle Orientation off";
        } else {
            changeProp(inode, 3, "off");
            myElem.setAttribute("class", "fakelinkoff");
            myElem.title = "Toggle Orientation on";
        }
    }

    myTree.update();
};


//
// callback when "CancelStepThru" is pressed in Tree
//
var cancelStepThru = function () {
    // alert("in cancelStepThru()");

    wv.curStep = 0;
    var button = document.getElementById("stepThruBtn");
    button["innerHTML"] = "StepThru";

    browserToServer("nextStep|0|");
};


//
// callback when "DisplayType" is pressed
//
var modifyDisplayType = function (e) {
    // alert("in modifyDisplayType(e="+e+")");

    if (checkIfWithBall()) {
        var ptype = prompt("Enter display type:\n"+
                           "   0  monochrome\n"+
                           "   1  normalized U parameter\n"+
                           "   2  normalized V parameter\n"+
                           "   3  minimum curvature\n"+
                           "   4  maximum curvature\n"+
                           "   5  Gaussian curvature\n"+
                           "   6  normals\n"+
                           "  10  Erep", "0");

        if (ptype === null) {
            return;
        } else if (isNaN(ptype)) {
            alert("Illegal number format, so no change");
            return;
        } else {
            wv.plotType = Number(ptype);

            if (wv.plotType > 0 && wv.plotType < 7) {
                setKeyLimits(null);
            } else {
                // send the limits back to the server
                browserToServer("setLims|"+wv.plotType+"|"+wv.loLimit+"|"+wv.upLimit+"|");
            }
        }
    }
};


//
// callback when "DisplayFilter" is pressed
//
var modifyDisplayFilter = function (e) {
    // alert("in modifyDisplayFilter(e="+e+")")

    if (checkIfFree() === false) {
        return;
    }

    if (typeof e != "object") {
        var attrName = e;
    } else {
        var attrName = prompt("Enter Attribute name (or * for all or ? for list)");
    }
    if (attrName === null) {
        displayFilterOff();
        return;
    } else if (attrName.length <= 0) {
        displayFilterOff();
        return;
    } else if (attrName == "?") {
        var attrNameList = "Attribute names found: ";

        for (var sgItem in wv.sgData) {
            var attrs = wv.sgData[sgItem];
            if (attrs !== undefined) {
                for (var i = 0; i < attrs.length; i+=2) {
                    if        (attrs[i].charAt(0) == "_") {
                        continue;
                    } else if (attrs[i].charAt(0) == ".") {
                        continue;
                    } else if (attrNameList.indexOf(" "+attrs[i]+",") < 0) {
                        attrNameList += attrs[i]+", ";
                    }
                }
            }
        }

        for (var gprim in wv.sceneGraph) {
            var attrs = wv.sgData[gprim];
            if (attrs !== undefined) {
                for (var i = 0; i < attrs.length; i+=2) {
                    if        (attrs[i].charAt(0) == "_") {
                        continue;
                    } else if (attrs[i].charAt(0) == ".") {
                        continue;
                    } else if (attrNameList.indexOf(" "+attrs[i]+",") < 0) {
                        attrNameList += attrs[i]+", ";
                    }
                }
            }
        }

        alert(attrNameList);
        modifyDisplayFilter(null);
        return;
    } else {
        var attrValue = prompt("Enter Attribute value (or * for all or ? for list)");
        if (attrValue === null) {
            displayFilterOff();
            return;
        } else if (attrValue.length <= 0) {
            displayFilterOff();
            return;
        } else if (attrValue == "?") {
            var attrValueList = "Attribute values for \"" + attrName + "\" found: ";

            for (var sgItem in wv.sgData) {
                var attrs = wv.sgData[sgItem];
                if (attrs !== undefined) {
                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrName == "*" || attrs[i] == attrName) {
                            if (attrValueList.indexOf(" "+attrs[i+1].trim()+",") < 0) {
                                attrValueList += attrs[i+1].trim()+", ";
                            }
                        }
                    }
                }
            }

            for (var gprim in wv.sceneGraph) {
                var attrs = wv.sgData[gprim];
                if (attrs !== undefined) {
                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrName == "*" || attrs[i] == attrName) {
                            if (attrValueList.indexOf(" "+attrs[i+1].trim()+",") < 0) {
                                attrValueList += attrs[i+1].trim()+", ";
                            }
                        }
                    }
                }
            }

            alert(attrValueList);
            modifyDisplayFilter(attrName);
            return;
        }

        if (attrName == "*" && attrValue == "*") {
            displayFilterOff();
            return;
        }

        // make all gprims transparent
        for (var gprim in wv.sceneGraph) {
            wv.sceneGraph[gprim].attrs |= wv.plotAttrs.TRANSPARENT;
        }

        // loop through all Bodys.  if it matches the filter,
        //    make all of its graphics primitives non-transparent
        for (var sgItem in wv.sgData) {
            var foo = sgItem.trim().split(" ");

            if ((foo.length == 1) || (foo.length == 2 && foo[0] == "Body")) {
                var attrs = wv.sgData[sgItem];
                for (var i = 0; i < attrs.length; i+=2) {
                    if        (attrs[i].trim() == attrName && attrs[i+1].trim() == attrValue) {
                        for (gprim in wv.sceneGraph) {
                            var bas = gprim.split(" ");
                            if ((foo.length == 1 && foo[0] == bas[0]) ||
                                (foo.length == 2 && foo[1] == bas[1])   ) {
                                wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            }
                        }
                    } else if ("*"             == attrName && attrs[i+1].trim() == attrValue) {
                        for (gprim in wv.sceneGraph) {
                            var bas = gprim.split(" ");
                            if ((foo.length == 1 && foo[0] == bas[0]) ||
                                (foo.length == 2 && foo[1] == bas[1])   ) {
                                wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            }
                        }
                    } else if (attrs[i].trim() == attrName && "*"               == attrValue) {
                        for (gprim in wv.sceneGraph) {
                            var bas = gprim.split(" ");
                            if ((foo.length == 1 && foo[0] == bas[0]) ||
                                (foo.length == 2 && foo[1] == bas[1])   ) {
                                wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                            }
                        }
                    }
                }
            }
        }


        // loop through all graphics primitives.  if it matches the filter,
        //    make it non-transparent
        if (isNaN(Number(attrValue)) == true) {
            for (gprim in wv.sceneGraph) {
                try {
                    var attrs = wv.sgData[gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        if        (attrs[i].trim() == attrName && attrs[i+1].trim() == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        } else if ("*"             == attrName && attrs[i+1].trim() == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        } else if (attrs[i].trim() == attrName && "*"               == attrValue) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        }
                    }
                } catch (x) {
                }
            }
            postMessage("Display filtered to \""+attrName+"\" \""+attrValue+"\"");
        } else {
            attrValue = Number(attrValue);
            var scale = 1.0 / Math.max(Math.abs(attrValue), 1e-6);

            for (gprim in wv.sceneGraph) {
                try {
                    var attrs = wv.sgData[gprim];
                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrs[i].trim() == attrName && scale*Math.abs(Number(attrs[i+1].trim())-attrValue) < 1e-3) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        } else if ("*"      == attrName && scale*Math.abs(Number(attrs[i+1].trim())-attrValue) < 1e-3) {
                            wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        }
                    }
                } catch (x) {
                }
            }
            postMessage("Display filtered to \""+attrName+"\" "+attrValue+"");
        }
    }

    myTree.update();
    wv.sceneUpd = 1;
};


//
// turn display filtering off
//
var displayFilterOff = function () {
    // alert("in displayFilterOff()");

    for (var gprim in wv.sceneGraph) {
        wv.sceneGraph[gprim].attrs &= ~wv.plotAttrs.TRANSPARENT;
    }

    postMessage("Display filtering off");

    myTree.update();
    wv.sceneUpd = 1;
};


//
// create wireframes from _wireframe attributes on Bodys
//
var createWireframes = function () {
    // alert("in createWireframes()");

    // set up wireframe data
    wv.wireframe = [];
    var iwire = 0;
    var names = wv.bodynames.split("|");
    for (var iname = 0; iname < names.length; iname++) {
        if (names[iname].length <= 0) continue;

        try {
            var attrs = wv.sgData[names[iname]];
            for (var iattr = 0; iattr < attrs.length; iattr+= 2) {
                if (attrs[iattr] == "_wireframe") {
                    var bbox = attrs[iattr+1].split(" ");
                    wv.wireframe[iwire++] = Number(bbox[1]);
                    wv.wireframe[iwire++] = Number(bbox[3]);
                    wv.wireframe[iwire++] = Number(bbox[5]);

                    wv.wireframe[iwire++] = Number(bbox[2]);
                    wv.wireframe[iwire++] = Number(bbox[3]);
                    wv.wireframe[iwire++] = Number(bbox[5]);

                    wv.wireframe[iwire++] = Number(bbox[1]);
                    wv.wireframe[iwire++] = Number(bbox[4]);
                    wv.wireframe[iwire++] = Number(bbox[5]);

                    wv.wireframe[iwire++] = Number(bbox[2]);
                    wv.wireframe[iwire++] = Number(bbox[4]);
                    wv.wireframe[iwire++] = Number(bbox[5]);

                    wv.wireframe[iwire++] = Number(bbox[1]);
                    wv.wireframe[iwire++] = Number(bbox[3]);
                    wv.wireframe[iwire++] = Number(bbox[6]);

                    wv.wireframe[iwire++] = Number(bbox[2]);
                    wv.wireframe[iwire++] = Number(bbox[3]);
                    wv.wireframe[iwire++] = Number(bbox[6]);

                    wv.wireframe[iwire++] = Number(bbox[1]);
                    wv.wireframe[iwire++] = Number(bbox[4]);
                    wv.wireframe[iwire++] = Number(bbox[6]);

                    wv.wireframe[iwire++] = Number(bbox[2]);
                    wv.wireframe[iwire++] = Number(bbox[4]);
                    wv.wireframe[iwire++] = Number(bbox[6]);
                }
            }
        } catch (x) {

        }
    }
    if (iwire == 0) {
        wv.wireframe = undefined;
    }
}


//
// function to draw wireframes (if they exist)
//
var drawWireframes = function (gl) {

    // return immediately if we have not set focus yet or if we do not have wireframes
    if (wv.focus.length == 0) {
        return;
    } else if (wv.wireframe === undefined) {
        return;
    }

    // remember the mvpMatrix (to be used when converting 3D to 2D coordinates)
    wv.wireMatrix.load(wv.mvpMatrix);

    // define the hexahedron associated with each wireframe
    var nout = wv.wireframe.length / 24;

    // set up the transformed coordinates for the wireframe
    var verts = new Float32Array(72*nout);
    for (var iout = 0; iout < nout; iout++) {

        // 0 - 1
        verts[72*iout   ] = (wv.wireframe[24*iout   ] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+ 1] = (wv.wireframe[24*iout+ 1] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+ 2] = (wv.wireframe[24*iout+ 2] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+ 3] = (wv.wireframe[24*iout+ 3] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+ 4] = (wv.wireframe[24*iout+ 4] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+ 5] = (wv.wireframe[24*iout+ 5] - wv.focus[2]) / wv.focus[3];

        // 2 - 3
        verts[72*iout+ 6] = (wv.wireframe[24*iout+ 6] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+ 7] = (wv.wireframe[24*iout+ 7] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+ 8] = (wv.wireframe[24*iout+ 8] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+ 9] = (wv.wireframe[24*iout+ 9] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+10] = (wv.wireframe[24*iout+10] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+11] = (wv.wireframe[24*iout+11] - wv.focus[2]) / wv.focus[3];

        // 4 - 5
        verts[72*iout+12] = (wv.wireframe[24*iout+12] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+13] = (wv.wireframe[24*iout+13] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+14] = (wv.wireframe[24*iout+14] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+15] = (wv.wireframe[24*iout+15] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+16] = (wv.wireframe[24*iout+16] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+17] = (wv.wireframe[24*iout+17] - wv.focus[2]) / wv.focus[3];

        // 6 - 7
        verts[72*iout+18] = (wv.wireframe[24*iout+18] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+19] = (wv.wireframe[24*iout+19] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+20] = (wv.wireframe[24*iout+20] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+21] = (wv.wireframe[24*iout+21] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+22] = (wv.wireframe[24*iout+22] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+23] = (wv.wireframe[24*iout+23] - wv.focus[2]) / wv.focus[3];

        // 0 - 2
        verts[72*iout+24] = (wv.wireframe[24*iout   ] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+25] = (wv.wireframe[24*iout+ 1] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+26] = (wv.wireframe[24*iout+ 2] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+27] = (wv.wireframe[24*iout+ 6] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+28] = (wv.wireframe[24*iout+ 7] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+29] = (wv.wireframe[24*iout+ 8] - wv.focus[2]) / wv.focus[3];

        // 1 - 3
        verts[72*iout+30] = (wv.wireframe[24*iout+ 3] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+31] = (wv.wireframe[24*iout+ 4] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+32] = (wv.wireframe[24*iout+ 5] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+33] = (wv.wireframe[24*iout+ 9] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+34] = (wv.wireframe[24*iout+10] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+35] = (wv.wireframe[24*iout+11] - wv.focus[2]) / wv.focus[3];

        // 4 - 6
        verts[72*iout+36] = (wv.wireframe[24*iout+12] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+37] = (wv.wireframe[24*iout+13] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+38] = (wv.wireframe[24*iout+14] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+39] = (wv.wireframe[24*iout+18] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+40] = (wv.wireframe[24*iout+19] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+41] = (wv.wireframe[24*iout+20] - wv.focus[2]) / wv.focus[3];

        // 5 - 7
        verts[72*iout+42] = (wv.wireframe[24*iout+15] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+43] = (wv.wireframe[24*iout+16] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+44] = (wv.wireframe[24*iout+17] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+45] = (wv.wireframe[24*iout+21] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+46] = (wv.wireframe[24*iout+22] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+47] = (wv.wireframe[24*iout+23] - wv.focus[2]) / wv.focus[3];

        // 0 - 4
        verts[72*iout+48] = (wv.wireframe[24*iout   ] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+49] = (wv.wireframe[24*iout+ 1] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+50] = (wv.wireframe[24*iout+ 2] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+51] = (wv.wireframe[24*iout+12] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+52] = (wv.wireframe[24*iout+13] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+53] = (wv.wireframe[24*iout+14] - wv.focus[2]) / wv.focus[3];

        // 1 - 5
        verts[72*iout+54] = (wv.wireframe[24*iout+ 3] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+55] = (wv.wireframe[24*iout+ 4] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+56] = (wv.wireframe[24*iout+ 5] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+57] = (wv.wireframe[24*iout+15] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+58] = (wv.wireframe[24*iout+16] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+59] = (wv.wireframe[24*iout+17] - wv.focus[2]) / wv.focus[3];

        // 2 - 6
        verts[72*iout+60] = (wv.wireframe[24*iout+ 6] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+61] = (wv.wireframe[24*iout+ 7] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+62] = (wv.wireframe[24*iout+ 8] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+63] = (wv.wireframe[24*iout+18] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+64] = (wv.wireframe[24*iout+19] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+65] = (wv.wireframe[24*iout+20] - wv.focus[2]) / wv.focus[3];

        // 3 - 7
        verts[72*iout+66] = (wv.wireframe[24*iout+ 9] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+67] = (wv.wireframe[24*iout+10] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+68] = (wv.wireframe[24*iout+11] - wv.focus[2]) / wv.focus[3];
        verts[72*iout+69] = (wv.wireframe[24*iout+21] - wv.focus[0]) / wv.focus[3];
        verts[72*iout+70] = (wv.wireframe[24*iout+22] - wv.focus[1]) / wv.focus[3];
        verts[72*iout+71] = (wv.wireframe[24*iout+23] - wv.focus[2]) / wv.focus[3];
    }

    // set up colors for the lines (orange)
    var colors = new Uint8Array(72*nout);
    for (iout = 0; iout < 24*nout; iout++) {
        colors[3*iout  ] = 255;
        colors[3*iout+1] =  97;
        colors[3*iout+2] =  71;
    }

    // draw the wireframes
    gl.disableVertexAttribArray(2);
    gl.uniform1f(wv.u_wLightLoc, 0.0);

    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, verts, gl.STATIC_DRAW);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(0);

    var cbuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
    gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
    gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.enableVertexAttribArray(1);

    gl.drawArrays(gl.LINES, 0, 24*nout);
    gl.deleteBuffer(buffer);
    gl.deleteBuffer(cbuf);
    gl.uniform1f(wv.u_wLightLoc, 1.0);
};


//
// function to draw lastPoint (if it has been defined)
//
var drawLastPoint = function (gl) {

    // return immediately if we have not set focus yet or if we do not have wireframes
    if (wv.focus.length == 0) {
        return;
    } else if (wv.lastPoint === undefined) {
        return;
    }

    // set up the transformed coordinates for the lastPoint
    var verts = new Float32Array(3);

    verts[0] = (wv.lastPoint[0] - wv.focus[0]) / wv.focus[3];
    verts[1] = (wv.lastPoint[1] - wv.focus[1]) / wv.focus[3];
    verts[2] = (wv.lastPoint[2] - wv.focus[2]) / wv.focus[3];

    // set up colors for the lines (orange)
    var colors = new Uint8Array(3);
    colors[0] = 255;
    colors[1] =  97;
    colors[2] =  71;

    // draw the wireframes
    gl.disableVertexAttribArray(2);
    gl.uniform1f(wv.u_wLightLoc, 0.0);
    gl.uniform1f(wv.u_pointSizeLoc, 5.0);

    var buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
    gl.bufferData(gl.ARRAY_BUFFER, verts, gl.STATIC_DRAW);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 0, 0);
    gl.enableVertexAttribArray(0);

    var cbuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, cbuf);
    gl.bufferData(gl.ARRAY_BUFFER, colors, gl.STATIC_DRAW);
    gl.vertexAttribPointer(1, 3, gl.UNSIGNED_BYTE, false, 0, 0);
    gl.enableVertexAttribArray(1);

    gl.drawArrays(gl.POINTS, 0, 1);
    gl.deleteBuffer(buffer);
    gl.deleteBuffer(cbuf);
    gl.uniform1f(wv.u_wLightLoc, 1.0);
};


//
// constructor for a Tree
//
var Tree = function (doc, treeId) {
    // alert("in Tree(doc="+doc+", treeId="+treeId+")");

    // remember the document
    this.document = doc;
    this.treeId   = treeId;

    // arrays to hold the Nodes
    this.name    = new Array();
    this.tooltip = new Array();
    this.gprim   = new Array();
    this.click   = new Array();
    this.parent  = new Array();
    this.child   = new Array();
    this.next    = new Array();
    this.nprop   = new Array();
    this.opened  = new Array();

    this.prop1  = new Array();
    this.cbck1  = new Array();
    this.prop2  = new Array();
    this.cbck2  = new Array();
    this.prop3  = new Array();
    this.cbck3  = new Array();

    // initialize Node=0 (the root)
    this.name[  0]  = "**root**";
    this.tooltip[0] = "";
    this.gprim[ 0]  = "";
    this.click[ 0]  = null;
    this.parent[0]  = -1;
    this.child[ 0]  = -1;
    this.next[  0]  = -1;
    this.nprop[ 0]  =  0;
    this.prop1[ 0]  = "";
    this.cbck1[ 0]  = null;
    this.prop2[ 0]  = "";
    this.cbck2[ 0]  = null;
    this.prop3[ 0]  = "";
    this.cbck3[ 0]  = null;
    this.opened[0]  = +1;

    // initially, delay is not on
    this.delay = 0;

    // add methods
    this.addNode  = TreeAddNode;
    this.expand   = TreeExpand;
    this.contract = TreeContract;
    this.prop     = TreeProp;
    this.clear    = TreeClear;
    this.build    = TreeBuild;
    this.update   = TreeUpdate;
};


//
// add a Node to the Tree
//
var TreeAddNode = function (iparent, name, tooltip, gprim, click,
                            prop1, cbck1,
                            prop2, cbck2,
                            prop3, cbck3) {
    // alert("in TreeAddNode(iparent="+iparent+", name="+name+", tooltip="+tooltip+", gprim="+gprim+
    //                       ", click="+click+
    //                       ", prop1="+prop1+", cbck1="+cbck1+
    //                       ", prop2="+prop2+", cbck2="+cbck2+
    //                       ", prop3="+prop3+", cbck3="+cbck3+")");

    // validate the input
    if (iparent < 0 || iparent >= this.name.length) {
        alert("iparent="+iparent+" is out of range");
        return;
    }

    // find the next Node index
    var inode = this.name.length;

    // store this Node's values
    this.name[   inode] = name;
    this.tooltip[inode] = tooltip;
    this.gprim[  inode] = gprim;
    this.click[  inode] = click;
    this.parent[ inode] = iparent;
    this.child[  inode] = -1;
    this.next[   inode] = -1;
    this.nprop[  inode] =  0;
    this.opened[ inode] =  0;

    // store the properties
    if (prop1 !== undefined) {
        this.nprop[inode] = 1;
        this.prop1[inode] = prop1;
        this.cbck1[inode] = cbck1;
    }

    if (prop2 !== undefined) {
        this.nprop[inode] = 2;
        this.prop2[inode] = prop2;
        this.cbck2[inode] = cbck2;
    }

    if (prop3 !== undefined) {
        this.nprop[inode] = 3;
        this.prop3[inode] = prop3;
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
};


//
// build the Tree (ie, create the html table from the Nodes)
//
var TreeBuild = function () {
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
        newTR.appendChild(newTDcol1);

        var newTexta = doc.createTextNode("");
        newTDcol1.appendChild(newTexta);

        // table data "node"+inode+"col2"
        var newTDcol2 = doc.createElement("TD");
        newTDcol2.setAttribute("id", "node"+inode+"col2");
        if (this.click[inode] != null) {
            newTDcol2.className = "fakelinkcmenu";
            if (this.tooltip[inode].length > 0) {
                newTDcol2.title = this.tooltip[inode];
            }
        }
        newTR.appendChild(newTDcol2);

        var newTextb = doc.createTextNode(this.name[inode]);
        newTDcol2.appendChild(newTextb);

        var name = this.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

        for (var jbrch = 0; jbrch < brch.length; jbrch++) {
            if (brch[jbrch].name == name) {
                if (brch[jbrch].ileft == -3) {
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

        // table data "node"+inode+"col4"
        if (this.nprop[inode] > 1) {
            var newTDcol4 = doc.createElement("TD");
            newTDcol4.setAttribute("id", "node"+inode+"col4");
            if (this.cbck2[inode] != "") {
                newTDcol4.className = "fakelinkoff";
            }
            newTR.appendChild(newTDcol4);

            if (this.nprop[inode] == 2) {
                newTDcol4.setAttribute("colspan", "2");
            }

            var newTextd = doc.createTextNode(this.prop2[inode]);
            newTDcol4.appendChild(newTextd);
        }

        // table data "node"+inode+"col5"
        if (this.nprop[inode] > 2) {
            var newTDcol5 = doc.createElement("TD");
            newTDcol5.setAttribute("id", "node"+inode+"col5");
            if (this.cbck3[inode] != "") {
                newTDcol5.className = "fakelinkoff";
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
};


//
// clear the Tree
//
var TreeClear = function () {
    // alert("in TreeClear()");

    // remove all but the first Node
    this.name.splice(   1);
    this.tooltip.splice(1);
    this.gprim.splice(  1);
    this.click.splice(  1);
    this.parent.splice( 1);
    this.child.splice(  1);
    this.next.splice(   1);
    this.nprop.splice(  1);
    this.opened.splice( 1);

    this.prop1.splice(1);
    this.cbck1.splice(1);
    this.prop2.splice(1);
    this.cbck2.splice(1);
    this.prop3.splice(1);
    this.cbck3.splice(1);

    // reset the root Node
    this.parent[0] = -1;
    this.child[ 0] = -1;
    this.next[  0] = -1;
};


//
// contract a Node in the Tree
//
var TreeContract = function (inode) {
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
};


//
// expand a Node in the Tree
//
var TreeExpand = function (inode) {
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
};


//
// change a property of a Node
//
var TreeProp = function (inode, iprop, onoff) {
    // alert("in TreeProp(inode="+inode+", iprop="+iprop+", onoff="+onoff+")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode="+inode+" is out of range");
        return;
    } else if (onoff != "on" && onoff != "off") {
        alert("onoff="+onoff+" is not 'on' or 'off'");
        return;
    }

    if (this != myTree) {
        alert("this="+this+"   myTree="+myTree);
    }

    var thisNode = "";

    // set the property for inode
    if        (iprop == 1 && this.prop1[inode] == "Viz") {
        thisNode = this.document.getElementById("node"+inode+"col3");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ON;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ON;
            }
        }
    } else if (iprop == 1) {
    } else if (iprop == 2 && this.prop2[inode] == "Grd") {
        thisNode = this.document.getElementById("node"+inode+"col4");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.LINES;
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.POINTS;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.LINES;
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.POINTS;
            }
        }
    } else if (iprop == 2) {
    } else if (iprop == 3 && this.prop3[inode] == "Trn") {
        thisNode = this.document.getElementById("node"+inode+"col5");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.TRANSPARENT;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.TRANSPARENT;
            }
        }
    } else if (iprop == 3 && this.prop3[inode] == "Ori") {
        thisNode = this.document.getElementById("node"+inode+"col5");

        if (this.gprim[inode] != "") {
            if (onoff == "on") {
                wv.sceneGraph[this.gprim[inode]].attrs |=  wv.plotAttrs.ORIENTATION;
            } else {
                wv.sceneGraph[this.gprim[inode]].attrs &= ~wv.plotAttrs.ORIENTATION;
            }
        }
    } else if (iprop ==3) {
    } else {
        alert("iprop="+iprop+" is not 1, 2, or 3");
        return;
    }

    // update fakelinks in TreeWindow (needed when .attrs do not exist)
    if (thisNode != "") {
        if (onoff == "on") {
            thisNode.setAttribute("class", "fakelinkon");
            thisNode.title = "Toggle Orientation off";
        } else {
            thisNode.setAttribute("class", "fakelinkoff");
            thisNode.title = "Toggle Orientation on";
        }
    }

    // set property for inode's children
    for (var jnode = inode+1; jnode < this.parent.length; jnode++) {
        if (this.parent[jnode] == inode) {
            this.prop(jnode, iprop, onoff);
        }
    }

    wv.sceneUpd = 1;
};


//
// update the Tree (after build/expension/contraction/property-set)
//
var TreeUpdate = function () {
    // alert("in TreeUpdate()");

    // delay will be turned on during the sync operation
    if (this.delay == 1) {
        return;
    }

    var doc = this.document;

    // traverse the Nodes using depth-first search
    for (var inode = 1; inode < this.opened.length; inode++) {
        var element = doc.getElementById("node"+inode);

        // unhide the row
        element.style.display = "table-row";

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

                myElem.className = "fakelinkexpand";
                myElem.firstChild.nodeValue = "+";
                myElem.title   = "Expand";
                myElem.onclick = function () {
                    var thisNode = this.id.substring(4);
                    thisNode     = thisNode.substring(0,thisNode.length-4);
                    This.expand(thisNode);
                };

            } else {
                var myElem = doc.getElementById("node"+inode+"col1");
                var This   = this;

                myElem.className = "fakelinkexpand";
                myElem.firstChild.nodeValue = "-";
                myElem.title   = "Collapse";
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
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Vizability on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Vizability off";
                    }
                }
            }
        }

        if (this.nprop[inode] >= 2) {
            var myElem = doc.getElementById("node"+inode+"col4");
            myElem.onclick = this.cbck2[inode];

            if (this.prop2[inode] == "Grd") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Grid on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Grid off";
                    }
                }
            }
        }

        if (this.nprop[inode] >= 3) {
            var myElem = doc.getElementById("node"+inode+"col5");
            myElem.onclick = this.cbck3[inode];

            if (this.prop3[inode] == "Trn") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Transparency on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Transparency off";
                    }
                }
            } else if (this.prop3[inode] == "Ori") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Orientation on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Orientation off";
                    }
                }
            }
        }
    }
};


//
// activate the buildButton
//
var activateBuildButton = function () {
    // alert("in activateBuildButton()");

    var button = document.getElementById("solveButton");

    button["innerHTML"] = "Press to Re-build";
    button.style.backgroundColor = "#3FFF3F";
};


//
// callback when "onresize" event occurs (called by ESP.html)
//
// resize the frames (with special handling to width of tlframe and height of brframe)
//
var resizeFrames = function () {
    // alert("in resizeFrames()");

    var scrollSize = 24;

    // get the size of the client (minus amount to account for scrollbars)
    var body = document.getElementById("mainBody");
    var bodyWidth  = body.clientWidth  - scrollSize;
    var bodyHeight = body.clientHeight - scrollSize;

    // get the elements associated with the frames and the canvas
    var topleft = document.getElementById("tlframe");
    var butnfrm = document.getElementById("butnfrm");
    var treefrm = document.getElementById("treefrm");
    var toprite = document.getElementById("trframe");
    var botleft = document.getElementById("blframe");
    var botrite = document.getElementById("brframe");
    var canvas  = document.getElementById(wv.canvasID);
    var sketch  = document.getElementById("sketcher");
    var gloves  = document.getElementById("gloves");

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
    canvas.style.width  = canvWidth+"px";
    canvas.width        = canvWidth;
    sketch.style.width  = canvWidth+"px";
    sketch.width        = canvWidth;
    gloves.style.width  = canvWidth+"px";
    gloves.width        = canvWidth;

    // compute and set the heights of the frames
    //    (do not make brframe larger than 200px)
    var botmHeight = Math.round(0.20 * bodyHeight);
    if (botmHeight > 200)   botmHeight = 200;
    var  topHeight = bodyHeight - botmHeight;
    var canvHeight =  topHeight - scrollSize - 5;
    var keyHeight  = botmHeight - 25;

    topleft.style.height =  topHeight+"px";
    treefrm.style.height = (topHeight-90)+"px";
    toprite.style.height =  topHeight+"px";
    botleft.style.height = botmHeight+"px";
    botrite.style.height = botmHeight+"px";
    canvas.style.height  = canvHeight+"px";
    canvas.height        = canvHeight;
    sketch.style.height  = canvHeight+"px";
    sketch.height        = canvHeight;
    gloves.style.height  = canvHeight+"px";
    gloves.height        = canvHeight;

    // set up csm editor window
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.setSize(null, topHeight-100);
    }

    // set up canvas associated with key
    if (wv.canvasKY !== undefined) {
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
};


//
// callback when "unload" event occurs (called by ESP.html)
//
// inform server that browser was closed
//
var browserClosing = function () {
    // alert("in browserClosing()");

    browserToServer("userName|"+wv.myName+"|*closed*|");

    // the following is needed to ensure that the browser informs
    //    the server that it is being closed before actually closing
    for (var i = 1000; i >= 0; i--) {
        postMessage("browser will die in: "+i);
    }
};


//
// change mode for trframe
//
var changeMode = function (newMode) {
    // alert("in changeMode(newMode="+newMode+")");

    var webViewer        = document.getElementById("WebViewer");
    var addBrchForm      = document.getElementById("addBrchForm");
    var editBrchForm     = document.getElementById("editBrchForm");
    var addBrchHeader    = document.getElementById("addBrchHeader");
    var editBrchHeader   = document.getElementById("editBrchHeader");
    var editValuForm     = document.getElementById("editValuForm");
    var showOutpmtrsForm = document.getElementById("showOutpmtrsForm");
    var addValuHeader    = document.getElementById("addValuHeader");
    var editValuHeader   = document.getElementById("editValuHeader");
    var editorForm       = document.getElementById("editorForm");
    var sketcherForm     = document.getElementById("sketcherForm");
    var glovesForm       = document.getElementById("glovesForm");

    var wvKey            = document.getElementById("WVkey");
    var sketcherStatus   = document.getElementById("sketcherStatus");
    var timStatus        = document.getElementById("timStatus");
    var ESPlogo          = document.getElementById("ESPlogo");

    if (newMode < 0) {
        wv.usingMain = 0;

        // used to cause buttons such as "cmdSave" not to be active
        wv.curMode = newMode;

        var button = document.getElementById("stepThruBtn");
        button["innerHTML"] = "StepThru";

        return;
    } else if (newMode == 0) {
        wv.usingMain = 1;

        webViewer.hidden        = false;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = true;
        ESPlogo.hidden          = false;

        wv.curMode   = 0;
        wv.curTool   = main;
        wv.curPmtr   = -1;
        wv.curBrch   = -1;
        wv.afterBrch = -1;
        wv.menuEvent = undefined;
        wv.keyPress  = -1;
    } else if (newMode == wv.curMode) {
        return;
    } else if (newMode == 1) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = false;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = true;
        ESPlogo.hidden          = false;

        // unselect all items
        var elements = document.getElementById("addBrchForm").elements;
        for (var ielem = 0; ielem < elements.length; ielem++) {
            elements[ielem].checked = false;
        }

        wv.curTool = main;
        wv.curMode = 1;
    } else if (newMode == 2) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = false;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = true;
        ESPlogo.hidden          = false;

        addBrchHeader.hidden    = false;
        editBrchHeader.hidden   = true;

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curTool = main;
        wv.curMode = 2;
    } else if (newMode == 3) {
        wv.singMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = false;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = true;
        ESPlogo.hidden          = false;

        addBrchHeader.hidden    = true;
        editBrchHeader.hidden   = false;

        showBrchArgs();

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curTool = main;
        wv.curMode = 3;
    } else if (newMode == 4) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = false;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = true;
        ESPlogo.hidden          = false;

        addValuHeader.hidden    = false;
        editValuHeader.hidden   = true;

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curTool = main;
        wv.curMode = 4;
    } else if (newMode == 5) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = false;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        ESPlogo.hidden          = false;

        addValuHeader.hidden    = true;
        timStatus.hidden        = true;
        editValuHeader.hidden   = false;

        if (wv.getFocus !== undefined) {
            wv.getFocus.focus();
            wv.getFocus.select();
            wv.getFocus = undefined;
        }

        wv.curTool = main;
        wv.curMode = 5;
    } else if (newMode == 6) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = false;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        ESPlogo.hidden          = false;

        wv.curTool = main;
        wv.curMode = 6;
    } else if (newMode == 7) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = false;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = true;
        ESPlogo.hidden          = false;

        wv.curTool = main;
        wv.curMode = 7;
    } else if (newMode == 8) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = false;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = false;
        timStatus.hidden        = true;
        ESPlogo.hidden          = true;

        wv.curTool = sketch;
        wv.curMode = 8;
    } else if (newMode == 9) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = false;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = false;
        ESPlogo.hidden          = true;

        wv.curTool = gloves;
        wv.curMode = 9;
    } else if (newMode == 10) {
        wv.usingMain = 1;

        webViewer.hidden        = false;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = false;
        ESPlogo.hidden          = true;

        wv.curTool = ereped;
        wv.curMode = 10;
    } else if (newMode == 11) {
        wv.usingMain = 1;

        webViewer.hidden        = false;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = false;
        ESPlogo.hidden          = true;

        wv.curTool = plugs;
        wv.curMode = 11;
    } else if (newMode == 12) {
        wv.usingMain = 1;

        webViewer.hidden        = false;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = false;
        ESPlogo.hidden          = true;

        wv.curTool = pyscript;
        wv.curMode = 12;
    } else if (newMode == 13) {
        wv.usingMain = 1;

        webViewer.hidden        = false;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = true;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = true;
        timStatus.hidden        = false;
        ESPlogo.hidden          = true;

        wv.curTool = mitten;
        wv.curMode = 13;
    } else if (newMode == 14) {
        wv.usingMain = 0;

        webViewer.hidden        = true;
        addBrchForm.hidden      = true;
        editBrchForm.hidden     = true;
        editValuForm.hidden     = true;
        showOutpmtrsForm.hidden = true;
        editorForm.hidden       = true;
        sketcherForm.hidden     = false;
        glovesForm.hidden       = true;
        wvKey.hidden            = true;
        sketcherStatus.hidden   = false;
        timStatus.hidden        = true;
        ESPlogo.hidden          = true;

        wv.curTool = plotter;
        wv.curMode = 14;
    } else {
        alert("Bad new mode = "+newMode);
    }

    if        (wv.curMode == 0) {
        document.getElementById("toolMenuBtn").hidden = false;
        document.getElementById("doneMenuBtn").hidden = true;
    } else if (wv.curMode >= 8) {
        document.getElementById("toolMenuBtn").hidden = true;
        document.getElementById("doneMenuBtn").hidden = false;
    }
};


//
// rebuild the Tree Window
//
var rebuildTreeWindow = function (x) {
    // alert("in rebuildTreeWindow("+x+")");

    // do not let this routine be executed more than once at a time
    if (wv.buildTree > 0) {
        return;
    } else {
        wv.buildTree = 1;
    }

    // if there was a previous Tree, keep track of whether or not
    //    the Parameters, Branches, and Display was open
    var cvalsOpen = 0;
    var pmtr1Open = 0;
    var pmtr2Open = 0;
    var brchsOpen = 0;

    for (var inode = 0; inode < myTree.opened.length; inode++) {
        if (myTree.name[inode] == "Caps Value") {
            cvalsOpen = myTree.opened[inode];
        } else if (myTree.name[inode] == "Design Parameters") {
            pmtr1Open = myTree.opened[inode];
        } else if (myTree.name[inode] == "Local Variables") {
            pmtr2Open = myTree.opened[inode];
        } else if (myTree.name[inode] == "Branches") {
            brchsOpen = myTree.opened[inode];
        } else if (myTree.name[inode] == "Display") {
            break;
        }
    }

    // clear previous Nodes from the Tree
    myTree.clear();

    wv.bodynames = "|";

    // put the group headers into the Tree
    if (wv.server == "serveCAPS") {
        myTree.addNode(0, "Caps Values",   "Add a Caps Value",  "", addCval,
                       "ExpandAll", allCvals);
        var capsValues = myTree.name.length - 1;
    }

    myTree.addNode(0, "Design Parameters", "Add a Parameter",   "", addPmtr,
                   "ExpandAll", allPmtrs);
    var designParameters = myTree.name.length - 1;

    myTree.addNode(0, "Local Variables",   "Show Outpmtrs",     "", showOutpmtrs,
                   "ExpandAll", allLocals);
    var localVariables = myTree.name.length - 1;

    if (wv.server != "serveCAPS") {
        myTree.addNode(0, "Branches",          "Add Branch at end", "", addBrch,
                       "ExpandAll", allBrchs);
        var branches = myTree.name.length - 1;
    }

    if (wv.curStep == 0) {
        myTree.addNode(0, "Display",       "Change display",    "", chgDisplay,
                       "Viz", toggleViz,
                       "Grd", toggleGrd);
    } else {
        myTree.addNode(0, "Display");
    }
    var display = myTree.name.length - 1;

    myTree.addNode(localVariables, "@-parameters",      "",     "", null   );
    var atParameters = myTree.name.length - 1;

    myTree.addNode(localVariables, "@@-parameters",     "",     "", null   );
    var atAtParameters = myTree.name.length - 1;

    // put the Caps Variables into the Tree
    if (wv.server == "serveCAPS") {
        for (var icval = 0; icval < cval.length; icval++) {
            var name   = "\u00a0\u00a0"+cval[icval].name;
            var nrow   =                cval[icval].nrow;
            var ncol   =                cval[icval].ncol;
            var value  =                cval[icval].value[0];
            var parent = capsValues;

            if (cval[icval].name.search(/^:.*/)   >= 0 ||
                cval[icval].name.search(/.*:$/)   >= 0 ||
                cval[icval].name.search(/.*::.*/) >= 0 ||
                cval[icval].name.indexOf(":")     <= 0   ) {
                parent = capsValues;
            } else {
                var parts = cval[icval].name.split(":");

                // make sure all prefixes are in Tree
                for (var iii = 0; iii < parts.length-1; iii++) {
                    var found = 0;
                    for (var jjj = 0; jjj < myTree.name.length; jjj++) {
                        if (myTree.name[jjj].replace(/\u00a0/g, "") == parts[iii]+":" &&
                            myTree.parent[jjj]                      == parent            ) {
                            parent = jjj;
                            found = 1;
                            break;
                        }
                    }
                    if (found == 0) {
                        var temp = parts[iii] + ":";
                        for (jjj = 0; jjj <= iii; jjj++) {
                            temp = "\u00a0\u00a0" + temp;
                        }
                        myTree.addNode(parent, temp, "", null);
                        parent = myTree.name.length - 1;
                    }
                }
                name = parts[parts.length-1];
                for (iii = 0; iii < parts.length; iii++) {
                    name = "\u00a0\u00a0" + name;
                }
            }

            if (nrow > 1 || ncol > 1) {
                value = "["+nrow+"x"+ncol+"]";
            }

            myTree.addNode(parent, name, "Edit Caps Value", "", editCval,
                           ""+value, "");
        }

        wv.cvalStat = -2;

        // open the Caps Values (if they were open before the Tree was rebuilt)
        if (cvalsOpen == 1) {
            myTree.opened[capsValues] = 1;
        }
    }

    // put the Design Parameters into the Tree
    for (var ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (pmtr[ipmtr].type == OCSM_DESPMTR ||
            pmtr[ipmtr].type == OCSM_CFGPMTR   ) {
            var name   = "\u00a0\u00a0"+pmtr[ipmtr].name;
            var nrow   =                pmtr[ipmtr].nrow;
            var ncol   =                pmtr[ipmtr].ncol;
            var value  =                pmtr[ipmtr].value[0];
            var parent = designParameters;

            if (pmtr[ipmtr].name.search(/^:.*/)   >= 0 ||
                pmtr[ipmtr].name.search(/.*:$/)   >= 0 ||
                pmtr[ipmtr].name.search(/.*::.*/) >= 0 ||
                pmtr[ipmtr].name.indexOf(":")     <= 0   ) {
                parent = designParameters;
            } else {
                var parts = pmtr[ipmtr].name.split(":");

                // make sure all prefixes are in Tree
                for (var iii = 0; iii < parts.length-1; iii++) {
                    var found = 0;
                    for (var jjj = 0; jjj < myTree.name.length; jjj++) {
                        if (myTree.name[jjj].replace(/\u00a0/g, "") == parts[iii]+":" &&
                            myTree.parent[jjj]                      == parent            ) {
                            parent = jjj;
                            found = 1;
                            break;
                        }
                    }
                    if (found == 0) {
                        var temp = parts[iii] + ":";
                        for (jjj = 0; jjj <= iii; jjj++) {
                            temp = "\u00a0\u00a0" + temp;
                        }
                        myTree.addNode(parent, temp, "", null);
                        parent = myTree.name.length - 1;
                    }
                }
                name = parts[parts.length-1];
                for (iii = 0; iii < parts.length; iii++) {
                    name = "\u00a0\u00a0" + name;
                }
            }

            if (nrow > 1 || ncol > 1) {
                value = "["+nrow+"x"+ncol+"]";
            }

            if (pmtr[ipmtr].dot[0] != 0) {
                name = "^" + name;
            }

            myTree.addNode(parent, name, "Edit Parameter", "", editPmtr,
                               ""+value, "");
        }
    }

    wv.pmtrStat = -2;

    // open the Design Parameters (if they were open before the Tree was rebuilt)
    if (pmtr1Open == 1) {
        myTree.opened[designParameters] = 1;
    }

    // put the Local Variables into the Tree
    for (ipmtr = 0; ipmtr < pmtr.length; ipmtr++) {
        if (pmtr[ipmtr].type != OCSM_DESPMTR &&
            pmtr[ipmtr].type != OCSM_CFGPMTR   ) {
            var name  = "\u00a0\u00a0"+pmtr[ipmtr].name;
            var nrow  =                pmtr[ipmtr].nrow;
            var ncol  =                pmtr[ipmtr].ncol;
            var value =                pmtr[ipmtr].value[0];

            if (nrow > 1 || ncol > 1) {
                value = "["+nrow+"x"+ncol+"]";
            }

            if (pmtr[ipmtr].name[0] == "@" && pmtr[ipmtr].name[1] == "@") {
                myTree.addNode(atAtParameters, name, "", "", null, ""+value, "");
            } else if (pmtr[ipmtr].name[0] == "@") {
                myTree.addNode(atParameters, name, "", "", null, ""+value, "");
            } else {
                myTree.addNode(localVariables, name, "", "", null, ""+value, "");
            }
            var inode = myTree.name.length - 1;

            var indx = 0;
            if (nrow > 1 || ncol > 1) {
                for (var irow = 0; irow < nrow; irow++) {
                    for (var icol = 0; icol < ncol; icol++) {
                        name  = "\u00a0\u00a0\u00a0\u00a0["+(irow+1)+","+(icol+1)+"]";
                        value = pmtr[ipmtr].value[indx++];

                        myTree.addNode(inode, name, "", "", null,
                                       "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0"+value, "");
                    }
                }
            }
        }
    }

    wv.pmtrStat = -2;

    // open the Local Variables (if they were open before the Tree was rebuilt)
    if (pmtr2Open == 1) {
        myTree.opened[localVariables] = 1;
    }

    // put the Branches into the Tree
    if (wv.server != "serveCAPS") {
        var parents = [branches];
        for (var ibrch = 0; ibrch < brch.length; ibrch++) {
            var name  = "\u00a0\u00a0";
            for (var indent = 0; indent < brch[ibrch].indent; indent++) {
                name = name+">";
            }
            name = name+brch[ibrch].name;

            if (ibrch == 0) {
            } else if (brch[ibrch].indent > brch[ibrch-1].indent) {
                parents.push(myTree.name.length-1);
            } else if (brch[ibrch].indent < brch[ibrch-1].indent) {
                parents.pop();
            }

            var type = brch[ibrch].type;
            var actv;
            if (ibrch == wv.curStep-1) {
                actv = "<<step<<";
            } else if (ibrch >= wv.builtTo) {
                actv = "skipped";
            } else if (brch[ibrch].actv == OCSM_SUPPRESSED) {
                actv = "suppressed";
            } else if (brch[ibrch].actv == OCSM_INACTIVE) {
                actv = "inactive";
            } else if (brch[ibrch].actv == OCSM_DEFERRED) {
                actv = "deferred";
            } else {
                actv = "";
            }

            myTree.addNode(parents[parents.length-1], name, "Edit/del/add-after Branch", "", editBrch,
                           type, "",
                           actv, "");
        }
        parents = undefined;

        wv.brchStat = -2;

        // open the Branches (if they were open before the Tree was rebuilt)
        if (brchsOpen == 1) {
            myTree.opened[branches] = 1;
        }
    }

    // put the Display attributes into the Tree
    var patchesNode = -1;    // tree Node that will contain the Patches
    for (var gprim in wv.sceneGraph) {
        if (wv.curStep > 0) {
            myTree.addNode(display, "\u00a0\u00a0CancelStepThru", "Cancel StepThru mode", null, cancelStepThru);
            break;
        }

        // parse the name
        var matches = gprim.split(" ");

        var iface = undefined;
        var iedge = undefined;
        var inode = undefined;
        var csys  = undefined;

        if        (matches[0] == "Axes") {
            myTree.addNode(display, "\u00a0\u00a0Axes", "", gprim, null,
                           "Viz", toggleViz);
            myTree.addNode(display, "\u00a0\u00a0DisplayType", "Modify display type",     null, modifyDisplayType);
            myTree.addNode(display, "\u00a0\u00a0DisplayFilter", "Modify display filter", null, modifyDisplayFilter);
            continue;

        // processing for a Patch: "Patch m @I=n"
        } else if (matches.length == 3 && matches[0] == "Patch" &&
                   (matches[2].includes("@I=") ||
                    matches[2].includes("@J=") ||
                    matches[2].includes("@K=")   )                 ) {
            if (patchesNode == -1) {
                myTree.addNode(display, "\u00a0\u00a0Patches", "", gprim, null,
                              "Viz", toggleViz);
                patchesNode = myTree.name.length - 1;
            }
            myTree.addNode(patchesNode, "\u00a0\u00a0\u00a0\u00a0"+gprim, "", gprim, null,
                           "Viz", toggleViz);
            continue;  // no further processing for this gprim

        // processing for plotdata: "PlotCP: body:face"
        } else if (matches[0] == "PlotCP:") {
            myTree.addNode(display, "\u00a0\u00a0PlotCP: "+matches[1], "", gprim, null,
                          "Viz", toggleViz);
            continue;  // no further processing for this gprim

        // processing for plotdata: "PlotPoints: name"
        } else if (matches[0] == "PlotPoints:") {
            myTree.addNode(display, "\u00a0\u00a0"+matches[1], "", gprim, null,
                          "Viz", toggleViz);
            continue;  // no further processing for this gprim

        // processing for plotdata: "PlotLine: name"
        } else if (matches[0] == "PlotLine:") {
            myTree.addNode(display, "\u00a0\u00a0"+matches[1], "", gprim, null,
                          "Viz", toggleViz);
            continue;  // no further processing for this gprim

        // processing for plotdata: "PlotTris: name"
        } else if (matches[0] == "PlotTris:") {
            myTree.addNode(display, "\u00a0\u00a0"+matches[1], "", gprim, null,
                          "Viz", toggleViz, "Grd", toggleGrd);
            continue;  // no further processing for this gprim

        // processing for plotdata: "PlotGrid: name"
        } else if (matches[0] == "PlotGrid:") {
            myTree.addNode(display, "\u00a0\u00a0"+matches[1], "", gprim, null,
                          "Viz", toggleViz);
            continue;  // no further processing for this gprim

        // processing when Body is explicitly named: "Bodyname"
        } else if (matches.length == 1) {
            var bodyName = matches[0];
//$$$            var iface    = -1;
//$$$            var iedge    = -1;
//$$$            var inode    = -1;
//$$$            var csys     = -1;

        // processing when Body is not explicitly named: "Body m"
        } else if (matches.length == 2 && matches[0] == "Body") {
            var bodyName = matches[0] + " " + matches[1];
//$$$            var iface    = -1;
//$$$            var iedge    = -1;
//$$$            var inode    = -1;
//$$$            var icsys    = -1;

        // processing when Body is explicitly named: "Bodyname Edge m"
        } else if (matches.length == 3) {
            var bodyName = matches[0];
            if        (matches[1] == "Face") {
                var iface = matches[2];
            } else if (matches[1] == "Edge") {
                var iedge = matches[2];
            } else if (matches[1] == "Node") {
                var inode = matches[2];
            } else if (matches[1] == "Csys") {
                var icsys = matches[2];
            }

        // processing when Body is not explicitly named: "Body m Edge n"
        } else if (matches.length == 4 && matches[0] == "Body") {
            var bodyName = matches[0] + " " + matches[1];
            if        (matches[2] == "Face") {
                var iface = matches[3];
            } else if (matches[2] == "Edge") {
                var iedge = matches[3];
            } else if (matches[2] == "Node") {
                var inode = matches[3];
            } else if (matches[2] == "Csys") {
                var icsys  = matches[3];
            }
        }

        // determine if Body does not exist by looking for a tree node
        //    whose name matches the Body's name and which is a descendant
        //    of the "Display" node
        var kbody = -1;
        for (var jnode = 1; jnode < myTree.name.length; jnode++) {
            if (myTree.name[jnode] == "\u00a0\u00a0"+bodyName) {
                var knode = jnode;
                while (knode != 0) {
                    if (knode == display) {             // "Display"
                        kbody = jnode;
                        break;
                    }
                    knode = myTree.parent[knode];
                }
            }
        }

        // if Body does not exist, create it and its Face, Edge, Node, and Csystem lists
        //    subnodes now
        var kface, kedge, knode, kcsys;
        if (kbody < 0) {
            myTree.addNode(display, "\u00a0\u00a0"+bodyName, "Show Body Attributes", "", showBodyAttrs,
                           "Viz", toggleViz,
                           "Grd", toggleGrd);
            kbody = myTree.name.length - 1;

            wv.bodynames += bodyName + "|";

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Faces", "", "", null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Trn", toggleTrn);
            kface = myTree.name.length - 1;

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Edges", "", "", null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Ori", toggleOri);
            kedge = myTree.name.length - 1;

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Nodes", "", "", null,
                           "Viz", toggleViz);
            knode = myTree.name.length - 1;

            myTree.addNode(kbody, "\u00a0\u00a0\u00a0\u00a0Csystems", "", "", null,
                           "Viz", toggleViz);
            kcsys = myTree.name.length - 1;

        // otherwise, get pointers to the face-group, edge-group, or node-group Nodes
        } else {
            kface = myTree.child[kbody];
            kedge = kface + 1;
            knode = kedge + 1;
            kcsys = knode + 1;
        }

        // make the Tree Node
        if        (iface !== undefined) {
            myTree.addNode(kface, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0face "+iface, "", gprim, null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Trn", toggleTrn);
        } else if (iedge !== undefined) {
            myTree.addNode(kedge, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0edge "+iedge, "", gprim, null,
                           "Viz", toggleViz,
                           "Grd", toggleGrd,
                           "Ori", toggleOri);
        } else if (inode !== undefined) {
            myTree.addNode(knode, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0node "+inode, "", gprim, null,
                           "Viz", toggleViz);
        } else if (icsys !== undefined) {
            myTree.addNode(kcsys, "\u00a0\u00a0\u00a0\u00a0\u00a0\u00a0"+icsys, "", gprim, null,
                           "Viz", toggleViz);
        }
    }

    // create wireframe data from _wireframe attributes
    createWireframes();

    // open the Display (by default)
    myTree.opened[display] = 1;

    // mark that we have (re-)built the Tree
    wv.sgUpdate = 0;

    // convert the abstract Tree Nodes into an HTML table
    myTree.build();

    wv.buildTree = 0;
};


//
// check if browser currently has the ball
//
var checkIfWithBall = function () {
    if (wv.myRole != 0) {
        if (confirm("This action is not permitted when you don't have the ball\n"+
                    "Do you want to request the ball?")) {
            browserToServer("message|"+wv.myName+" says: please pass me the ball|");
            alert("Try again once you have the ball");
        }
        return false;      // we do not have the ball and should try once we get it
    }

    return true;
};


//
// check if browser is not sync'd to browser with ball
//
var checkIfFree = function () {
    if (wv.myRole == 1) {
        if (confirm("This action is not permitted when you are sync'd\n"+
                    "Do you want to un-sync?") !== true) {
            return false;     // we are sync's and should not perform the action
        } else {
            wv.myRole = 2;

            var button = document.getElementById("collabMenuBtn");
            button.style.backgroundColor = null;
        }
    }

    return true;
};


//
// post a message into the brframe
//
var postMessage = function (mesg) {
    // alert("in postMessage(mesg="+mesg+")");

    if (wv.debugUI) {
        console.log("postMessage: "+mesg.substring(0,40));
    }

    var botm = document.getElementById("brframe");

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);
    botm.insertBefore(pre, botm.lastChild);

    pre.scrollIntoView();

    var text2 = botm["innerHTML"].replace(/\n<\/pre><pre>/g, "");
    botm["innerHTML"] = text2;
};


//
// load info into editBrchForm
//
var setupEditBrchForm = function () {
    // alert("in setupEditBrchForm()");

    var ibrch = wv.curBrch;
    var name  = brch[ibrch].name;
    var type  = brch[ibrch].type;
    var level = brch[ibrch].level;

    // return if within a UDC
    if (level > 0) {
        return 1;
    }

    var editBrchForm = document.getElementById("editBrchForm");

    // turn all arguments off (by default)
    document.getElementById("argName1").parentNode.style.display = "none";
    document.getElementById("argName2").parentNode.style.display = "none";
    document.getElementById("argName3").parentNode.style.display = "none";
    document.getElementById("argName4").parentNode.style.display = "none";
    document.getElementById("argName5").parentNode.style.display = "none";
    document.getElementById("argName6").parentNode.style.display = "none";
    document.getElementById("argName7").parentNode.style.display = "none";
    document.getElementById("argName8").parentNode.style.display = "none";
    document.getElementById("argName9").parentNode.style.display = "none";

    // turn all attributes/csystems off (by default)
    document.getElementById("attrName1").parentNode.style.display = "none";
    document.getElementById("attrName2").parentNode.style.display = "none";
    document.getElementById("attrName3").parentNode.style.display = "none";
    document.getElementById("attrName4").parentNode.style.display = "none";
    document.getElementById("attrName5").parentNode.style.display = "none";
    document.getElementById("attrName6").parentNode.style.display = "none";
    document.getElementById("attrName7").parentNode.style.display = "none";
    document.getElementById("attrName8").parentNode.style.display = "none";
    document.getElementById("attrName9").parentNode.style.display = "none";

    // start by looking at arguments
    document.getElementById("editBrchArgs" ).hidden = false;
    document.getElementById("editBrchAttrs").hidden = true;

    // fill in the Branch type and name
    document.getElementById("brchType").firstChild["data"] = type;
    editBrchForm.            brchName.value                = name;

    // by default, the "Enter Sketcher" button is turned off
    document.getElementById("EnterSketcher").style.display = "none";

    // set up arguments based upon the type
    var argList;
    var defValue;
    var suppress = 0;   // set to 1 if type can be suppressed

    if        (type == "applycsys") {
        argList  = ["$csysName", "ibody"];
        defValue = ["",          "0"    ];
    } else if (type == "arc") {
        argList  = ["xend", "yend", "zend", "dist", "$plane"];
        defValue = ["",     "",     "",     "",     "xy"    ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "assert") {
        argList  = ["arg1", "arg2", "toler", "verify"];
        defValue = ["",     "",     "0",     "0"     ];
    } else if (type == "bezier") {
        argList  = ["x", "y", "z"];
        defValue = ["",  "",  "" ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "blend") {
        argList  = ["begList", "endList", "reorder", "oneFace", "periodic"];
        defValue = ["0",       "0",       "0",       "0"      , "0"       ];
        suppress = 1;
    } else if (type == "box") {
        argList  = ["xmin", "ymin", "zmin", "dx", "dy", "dz"];
        defValue = ["",     "",     "",     "",   "",   ""  ];
        suppress = 1;
    } else if (type == "catbeg") {
        argList  = ["errCode"];
        defValue = [""       ];
    } else if (type == "catend") {
        argList  = [];
        defValue = [];
    } else if (type == "chamfer") {
        argList  = ["radius", "edgeList", "listStyle"];
        defValue = ["",       "0"       , "0"        ];
        suppress = 1;
    } else if (type == "cirarc") {
        argList  = ["xon", "yon", "zon", "xend", "yend", "zend"];
        defValue = ["",    "",    "",    "",     "",     ""    ];
    } else if (type == "combine") {
        argList  = ["toler"];
        defValue = ["0"    ];
    } else if (type == "cone") {
        argList  = ["xvrtx", "yvrtx", "zvrtx", "xbase", "ybase", "zbase", "radius"];
        defValue = ["",      "",      "",      "",      "",      "",      ""      ];
        suppress = 1;
    } else if (type == "connect") {
        argList  = ["faceList1", "faceList2", "edgelist1", "edgelist2", "toler"];
        defValue = ["",          ""         , "0",         "0",         "0"    ];
    } else if (type == "cylinder") {
        argList  = ["xbeg", "ybeg", "zbeg", "xend", "yend", "zend", "radius"];
        defValue = ["",     "",     "",     "",     "",     "",     ""      ];
        suppress = 1;
    } else if (type == "dimension") {
        argList  = ["$pmtrName", "nrow", "ncol"];
        defValue = ["",          "",     ""    ];
    } else if (type == "dump") {
        argList  = ["$filename", "remove", "tomark"];
        defValue = ["",          "0",      "0"     ];
    } else if (type == "else") {
        argList  = [];
        defValue = [];
    } else if (type == "elseif") {
        argList  = ["val1", "$op1", "val2", "$op2", "val3", "$op3", "val4"];
        defValue = ["",     "",     "",     "and",  "0",    "eq",   "0"   ];
//  } else if (type == "end") {
    } else if (type == "endif") {
        argList  = [];
        defValue = [];
    } else if (type == "evaluate") {
        argList  = ["$type", "arg1", "arg2", "arg3", "arg4", "arg5"];
        defValue = ["",      "",     "",     "",     "",     ""    ];
    } else if (type == "extract") {
        argList  = ["index"];
        defValue = [""     ];
    } else if (type == "extrude") {
        argList  = ["dx", "dy", "dz"];
        defValue = ["",   "",   ""  ];
        suppress = 1;
    } else if (type == "fillet") {
        argList  = ["radius", "edgeList", "listStyle"];
        defValue = ["",       "0"       , "0"        ];
        suppress = 1;
    } else if (type == "getattr") {
        argList  = ["$pmtrName", "attrID"];
        defValue = ["",          ""      ];
    } else if (type == "group") {
        argList  = [];
        defValue = [];
    } else if (type == "hollow") {
        argList  = ["thick",  "entList", "listStyle"];
        defValue = ["0",      "0"      , "0"        ];
        suppress = 1;
    } else if (type == "ifthen") {
        argList  = ["val1", "$op1", "val2", "$op2", "val3", "$op3", "val4"];
        defValue = ["",     "",     "",     "and",  "0",    "eq",   "0"   ];
    } else if (type == "import") {
        argList  = ["$filename", "bodynumber"];
        defValue = ["",          "1"         ];
        suppress = 1;
//  } else if (type == "interface") {
    } else if (type == "intersect") {
        argList  = ["$order", "index", "maxtol"];
        defValue = ["none",   "1",     "0"     ];
    } else if (type == "join") {
        argList  = ["toler"];
        defValue = ["0"    ];
    } else if (type == "linseg") {
        argList  = ["x", "y", "z"];
        defValue = ["",  "",  "" ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "loft") {
        alert("Consider using 'rule' or 'blend'");
        argList  = ["smooth"];
        defValue = [""      ];
        suppress = 1;
    } else if (type == "macbeg") {
        alert("'macbeg' is DEPRECATED");
        argList  = ["imacro"];
        defValue = [""      ];
    } else if (type == "macend") {
        alert("'macend' is DEPRECATED");
        argList  = [];
        defValue = [];
    } else if (type == "mark") {
        argList  = [];
        defValue = [];
    } else if (type == "message") {
        argList  = ["$text", "$schar"];
        defValue = ["",      "_"     ];
    } else if (type == "mirror") {
        argList  = ["nx", "ny", "nz", "dist"];
        defValue = ["",   "",   "",   "0"   ];
        suppress = 1;
    } else if (type == "patbeg") {
        argList  = ["$pmtrName", "ncopy"];
        defValue = ["",          ""     ];
    } else if (type == "patbreak") {
        argList  = ["expr"];
        defValue = [""    ];
    } else if (type == "patend") {
        argList  = [];
        defValue = [];
    } else if (type == "point") {
        argList  = ["xloc", "yloc", "zloc"];
        defValue = ["",     "",     ""    ];
        suppress = 1;
    } else if (type == "project") {
        argList  = ["x", "y", "z", "dx", "dy", "dz", "useEdges"];
        defValue = ["",  "",  "",  "",   "",   "",   "0"       ];
    } else if (type == "recall") {
        aleft("'recall' is DEPRECATED");
        argList  = ["imicro"];
        defValue = [""      ];
    } else if (type == "reorder") {
        argList  = ["ishift", "iflip"];
        defValue = ["",       "0"    ];
        suppress = 1;
    } else if (type == "restore") {
        argList  = ["$name", "index"];
        defValue = ["",      "0"    ];
    } else if (type == "revolve") {
        argList  = ["xorig", "yorig", "zorig", "dxaxis", "dyaxis", "dzaxis", "angDeg"];
        defValue = ["",      "",      "",      "",       "",       "",       ""      ];
        suppress = 1;
    } else if (type == "rotatex") {
        argList  = ["angDeg", "yaxis", "zaxis"];
        defValue = ["",       "0",     "0"    ];
        suppress = 1;
    } else if (type == "rotatey") {
        argList  = ["angDeg", "zaxis", "xaxis"];
        defValue = ["",       "0",     "0"    ];
        suppress = 1;
    } else if (type == "rotatez") {
        argList  = ["angDeg", "xaxis", "yaxis"];
        defValue = ["",       "0",     "0"     ];
        suppress = 1;
    } else if (type == "rule") {
        argList  = ["reorder", "periodic"];
        defValue = ["0"      , "0"       ];
        suppress = 1;
    } else if (type == "scale") {
        argList  = ["fact"];
        defValue = [""    ];
        suppress = 1;
    } else if (type == "select") {
        argList  = ["$type", "arg1", "arg2", "arg3", "arg4", "arg5", "arg6", "arg7", "arg8"];
        defValue = ["",      "",     "",     "",     "",     "",     "",     "",     ""    ];
    } else if (type == "set") {
        argList  = ["$pmtrName", "exprs"];
        defValue = ["",          ""     ];
    } else if (type == "skbeg") {
        argList  = ["x", "y", "z", "relative"];
        defValue = ["",  "",  "",  "1"       ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "skcon") {
        argList  = ["$type", "index1", "index2", "$value"];
        defValue = ["",      "",       "-1",     "0"     ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "skend") {
        argList  = ["wireOnly"];
        defValue = ["0"];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "skvar") {
        argList  = ["$type", "valList"];
        defValue = ["",      ""       ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "solbeg") {
        argList  = ["$varList"];
        defValue = [""        ];
    } else if (type == "solcon") {
        argList  = ["$expr"];
        defValue = [""     ];
    } else if (type == "solend") {
        argList  = [];
        defValue = [];
    } else if (type == "sphere") {
        argList  = ["xcent", "ycent", "zcent", "radius"];
        defValue = ["",      "",      "",      ""      ];
        suppress = 1;
    } else if (type == "spline") {
        argList  = ["x", "y", "z"];
        defValue = ["",  "",  "" ];
        document.getElementById("EnterSketcher").style.display = "inline";
    } else if (type == "sslope") {
        argList  = ["dx", "dy", "dz"];
        defValue = ["",  "",  "" ];
    } else if (type == "store") {
        argList  = ["$name", "index", "keep"];
        defValue = ["",      "0" ,    "0"   ];
    } else if (type == "subtract") {
        argList  = ["$order", "index", "maxtol"];
        defValue = ["none",   "1",     "0"     ];
    } else if (type == "sweep") {
        argList  = [];
        defValue = [];
        suppress = 1;
    } else if (type == "throw") {
        argList  = ["errCode"];
        defValue = [""       ];
    } else if (type == "torus") {
        argList  = ["xcent", "ycent", "zcent", "dxaxis", "dyaxis", "dzaxis", "majorRad", "minorRad"];
        defValue = ["",      "",      "",      "",       "",       "",       "",         ""        ];
        suppress = 1;
    } else if (type == "translate") {
        argList  = ["dx", "dy", "dz"];
        defValue = ["",   "",   ""  ];
        suppress = 1;
    } else if (type == "udparg") {
        argList  = ["$primtype", "$argName1", "argValue1", "$argName2", "argValue2", "$argName3", "argValue3", "$argName4", "argValue4"];
        defValue = ["",          "",          "",          "",          "",          "",          "",          "",          ""         ];
    } else if (type == "udprim") {
        argList  = ["$primtype", "$argName1", "argValue1", "$argName2", "argValue2", "$argName3", "argValue3", "$argName4", "argValue4"];
        defValue = ["",          "",          "",          "",          "",          "",          "",          "",          ""         ];
        suppress = 1;
    } else if (type == "union") {
        argList  = ["tomark", "trimList", "maxtol"];
        defValue = ["0",      "0",        "0"     ];
    } else {
        return 1;
    }

    // set up the activity
    if (suppress) {
        if        (brch[ibrch].actv == OCSM_ACTIVE) {     // Active
            editBrchForm.activity.value = OCSM_ACTIVE;
        } else if (brch[ibrch].actv == OCSM_SUPPRESSED) {     // Suppressed
            editBrchForm.activity.value = OCSM_SUPPRESSED;
        } else {                                  // Inactive or Deferred
            editBrchForm.activity.value = OCSM_ACTIVE;
        }
        editBrchForm.activity.style.display = "table-row";
    } else {
        editBrchForm.activity.style.display = "none";
    }

    // set up the number of Attributes
    if (wv.curMode != 1) {
        document.getElementById("numArgs").firstChild["data"]  = brch[ibrch].attrs.length;
    }

    // set up the arguments
    wv.numArgs = argList.length;

    if (wv.curMode != 1 && (type == "evaluate" || type == "select" || type == "udparg" || type == "udprim")) {
        if        (brch[ibrch].args[0] == "" || brch[ibrch].args[0] == "$") {
            wv.numArgs = 0;
        } else if (brch[ibrch].args[1] == "" || brch[ibrch].args[1] == "$") {
            wv.numArgs = 1;
        } else if (brch[ibrch].args[2] == "" || brch[ibrch].args[2] == "$") {
            wv.numArgs = 2;
        } else if (brch[ibrch].args[3] == "" || brch[ibrch].args[3] == "$") {
            wv.numArgs = 3;
        } else if (brch[ibrch].args[4] == "" || brch[ibrch].args[4] == "$") {
            wv.numArgs = 4;
        } else if (brch[ibrch].args[5] == "" || brch[ibrch].args[5] == "$") {
            wv.numArgs = 5;
        } else if (brch[ibrch].args[6] == "" || brch[ibrch].args[6] == "$") {
            wv.numArgs = 6;
        } else if (brch[ibrch].args[7] == "" || brch[ibrch].args[7] == "$") {
            wv.numArgs = 7;
        } else if (brch[ibrch].args[8] == "" || brch[ibrch].args[8] == "$") {
            wv.numArgs = 8;
        } else if (brch[ibrch].args[0] == "" || brch[ibrch].args[0] == "$") {
            wv.numArgs = 9;
        }
    }

    if (wv.numArgs >= 1) {
        document.getElementById("argName1").parentNode.style.display = "table-row";
        document.getElementById("argName1").firstChild["data"]       = argList[0];
        if (wv.curMode == 1 || brch[ibrch].args[0] === undefined) {
            editBrchForm.argValu1.value = defValue[0];
        } else if (argList[0].charAt(0) == "$") {
            editBrchForm.argValu1.value = brch[ibrch].args[0].substr(1,brch[ibrch].args[0].length);
        } else {
            editBrchForm.argValu1.value = brch[ibrch].args[0];
        }

        wv.getFocus = editBrchForm.argValu1;
    }
    if (wv.numArgs >= 2) {
        document.getElementById("argName2").parentNode.style.display = "table-row";
        document.getElementById("argName2").firstChild["data"]       = argList[1];
        if (wv.curMode == 1 || brch[ibrch].args[1] === undefined) {
            editBrchForm.argValu2.value = defValue[1];
        } else if (argList[1].charAt(0) == "$") {
            editBrchForm.argValu2.value = brch[ibrch].args[1].substr(1,brch[ibrch].args[1].length);
        } else {
            editBrchForm.argValu2.value = brch[ibrch].args[1];
        }
    }
    if (wv.numArgs >= 3) {
        document.getElementById("argName3").parentNode.style.display = "table-row";
        document.getElementById("argName3").firstChild["data"]       = argList[2];
        if (wv.curMode == 1 || brch[ibrch].args[2] === undefined) {
            editBrchForm.argValu3.value = defValue[2];
        } else if (argList[2].charAt(0) == "$") {
            editBrchForm.argValu3.value = brch[ibrch].args[2].substr(1,brch[ibrch].args[2].length);
        } else {
            editBrchForm.argValu3.value = brch[ibrch].args[2];
        }
    }
    if (wv.numArgs >= 4) {
        document.getElementById("argName4").parentNode.style.display = "table-row";
        document.getElementById("argName4").firstChild["data"]       = argList[3];
        if (wv.curMode == 1 || brch[ibrch].args[3] === undefined) {
            editBrchForm.argValu4.value = defValue[3];
        } else if (argList[3].charAt(0) == "$") {
            editBrchForm.argValu4.value = brch[ibrch].args[3].substr(1,brch[ibrch].args[3].length);
        } else {
            editBrchForm.argValu4.value = brch[ibrch].args[3];
        }
    }
    if (wv.numArgs >= 5) {
        document.getElementById("argName5").parentNode.style.display = "table-row";
        document.getElementById("argName5").firstChild["data"]       = argList[4];
        if (wv.curMode == 1 || brch[ibrch].args[4] === undefined) {
            editBrchForm.argValu5.value = defValue[4];
        } else if (argList[4].charAt(0) == "$") {
            editBrchForm.argValu5.value = brch[ibrch].args[4].substr(1,brch[ibrch].args[4].length);
        } else {
            editBrchForm.argValu5.value = brch[ibrch].args[4];
        }
    }
    if (wv.numArgs >= 6) {
        document.getElementById("argName6").parentNode.style.display = "table-row";
        document.getElementById("argName6").firstChild["data"]       = argList[5];
        if (wv.curMode == 1 || brch[ibrch].args[5] === undefined) {
            editBrchForm.argValu6.value = defValue[5];
        } else if (argList[5].charAt(0) == "$") {
            editBrchForm.argValu6.value = brch[ibrch].args[5].substr(1,brch[ibrch].args[5].length);
        } else {
            editBrchForm.argValu6.value = brch[ibrch].args[5];
        }
    }
    if (wv.numArgs >= 7) {
        document.getElementById("argName7").parentNode.style.display = "table-row";
        document.getElementById("argName7").firstChild["data"]       = argList[6];
        if (wv.curMode == 1 || brch[ibrch].args[6] === undefined) {
            editBrchForm.argValu7.value = defValue[6];
        } else if (argList[6].charAt(0) == "$") {
            editBrchForm.argValu7.value = brch[ibrch].args[6].substr(1,brch[ibrch].args[6].length);
        } else {
            editBrchForm.argValu7.value = brch[ibrch].args[6];
        }
    }
    if (wv.numArgs >= 8) {
        document.getElementById("argName8").parentNode.style.display = "table-row";
        document.getElementById("argName8").firstChild["data"]       = argList[7];
        if (wv.curMode == 1 || brch[ibrch].args[7] === undefined) {
            editBrchForm.argValu8.value = defValue[7];
        } else if (argList[7].charAt(0) == "$") {
            editBrchForm.argValu8.value = brch[ibrch].args[7].substr(1,brch[ibrch].args[7].length);
        } else {
            editBrchForm.argValu8.value = brch[ibrch].args[7];
        }
    }
    if (wv.numArgs >= 9) {
        document.getElementById("argName9").parentNode.style.display = "table-row";
        document.getElementById("argName9").firstChild["data"]       = argList[8];
        if (wv.curMode == 1 || brch[ibrch].args[8] === undefined) {
            editBrchForm.argValu9.value = defValue[8];
        } else if (argList[8].charAt(0) == "$") {
            editBrchForm.argValu9.value = brch[ibrch].args[8].substr(1,brch[ibrch].args[8].length);
        } else {
            editBrchForm.argValu9.value = brch[ibrch].args[8];
        }
    }

    if (wv.curMode != 1) {
        if (brch[ibrch].attrs.length > 0) {
            document.getElementById("attrName1").parentNode.style.display = "table-row";
            document.getElementById("attrName1").firstChild["data"]       = brch[ibrch].attrs[0][0];
            document.getElementById("attrType1").firstChild["data"]       = brch[ibrch].attrs[0][1];
            editBrchForm.attrValu1.value                                  = brch[ibrch].attrs[0][2];
        }
        if (brch[ibrch].attrs.length > 1) {
            document.getElementById("attrName2").parentNode.style.display = "table-row";
            document.getElementById("attrName2").firstChild["data"]       = brch[ibrch].attrs[1][0];
            document.getElementById("attrType2").firstChild["data"]       = brch[ibrch].attrs[1][1];
            editBrchForm.attrValu2.value                                  = brch[ibrch].attrs[1][2];
        }
        if (brch[ibrch].attrs.length > 2) {
            document.getElementById("attrName3").parentNode.style.display = "table-row";
            document.getElementById("attrName3").firstChild["data"]       = brch[ibrch].attrs[2][0];
            document.getElementById("attrType3").firstChild["data"]       = brch[ibrch].attrs[2][1];
            editBrchForm.attrValu3.value                                  = brch[ibrch].attrs[2][2];
        }
        if (brch[ibrch].attrs.length > 3) {
            document.getElementById("attrName4").parentNode.style.display = "table-row";
            document.getElementById("attrName4").firstChild["data"]       = brch[ibrch].attrs[3][0];
            document.getElementById("attrType4").firstChild["data"]       = brch[ibrch].attrs[3][1];
            editBrchForm.attrValu4.value                                  = brch[ibrch].attrs[3][2];
        }
        if (brch[ibrch].attrs.length > 4) {
            document.getElementById("attrName5").parentNode.style.display = "table-row";
            document.getElementById("attrName5").firstChild["data"]       = brch[ibrch].attrs[4][0];
            document.getElementById("attrType5").firstChild["data"]       = brch[ibrch].attrs[4][1];
            editBrchForm.attrValu5.value                                  = brch[ibrch].attrs[4][2];
        }
        if (brch[ibrch].attrs.length > 5) {
            document.getElementById("attrName6").parentNode.style.display = "table-row";
            document.getElementById("attrName6").firstChild["data"]       = brch[ibrch].attrs[5][0];
            document.getElementById("attrType6").firstChild["data"]       = brch[ibrch].attrs[5][1];
            editBrchForm.attrValu6.value                                  = brch[ibrch].attrs[5][2];
        }
        if (brch[ibrch].attrs.length > 6) {
            document.getElementById("attrName7").parentNode.style.display = "table-row";
            document.getElementById("attrName7").firstChild["data"]       = brch[ibrch].attrs[6][0];
            document.getElementById("attrType7").firstChild["data"]       = brch[ibrch].attrs[6][1];
            editBrchForm.attrValu7.value                                  = brch[ibrch].attrs[6][2];
        }
        if (brch[ibrch].attrs.length > 7) {
            document.getElementById("attrName8").parentNode.style.display = "table-row";
            document.getElementById("attrName8").firstChild["data"]       = brch[ibrch].attrs[7][0];
            document.getElementById("attrType8").firstChild["data"]       = brch[ibrch].attrs[7][1];
            editBrchForm.attrValu8.value                                  = brch[ibrch].attrs[7][2];
        }
        if (brch[ibrch].attrs.length > 8) {
            document.getElementById("attrName9").parentNode.style.display = "table-row";
            document.getElementById("attrName9").firstChild["data"]       = brch[ibrch].attrs[8][0];
            document.getElementById("attrType9").firstChild["data"]       = brch[ibrch].attrs[8][1];
            editBrchForm.attrValu9.value                                  = brch[ibrch].attrs[8][2];
        }
    }

    return 0;
};


//
// load info into editValuForm
//
var setupEditValuForm = function () {
    // alert("in setupEditValuForm()");

    // setting up for a Parameter
    if (wv.curPmtr >= 0) {
        var ipmtr = wv.curPmtr;
        var name  = pmtr[ipmtr].name;
        var nrow  = pmtr[ipmtr].nrow;
        var ncol  = pmtr[ipmtr].ncol;

        var editValuForm = document.getElementById("editValuForm");

        // fill in the Parameter name
        document.getElementById("valuName").firstChild["data"] = name;

        var valuValTable = document.getElementById("editValuValTable");

        // remove old table entries
        if (valuValTable) {
            var child1 = valuValTable.lastChild;
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
                valuValTable.removeChild(child1);
                child1 = valuValTable.lastChild;
            }
        }

        // build the table that will contain values
        for (var irow = 0; irow <= nrow; irow++) {
            var newTR = document.createElement("TR");
            valuValTable.appendChild(newTR);

            // fill the row
            if (irow == 0) {
                var newTD = document.createElement("TD");
                newTR.appendChild(newTD);

                var newText = document.createTextNode("");
                newTD.appendChild(newText);

                for (var icol = 1; icol <= ncol; icol++) {
                    newTD = document.createElement("TD");
                    newTR.appendChild(newTD);

                    newText = document.createTextNode("col\u00a0"+icol);
                    newTD.appendChild(newText);
                }
            } else{
                var newTD = document.createElement("TD");
                newTR.appendChild(newTD);

                var newText = document.createTextNode("row\u00a0"+irow);
                newTD.appendChild(newText);

                for (var icol = 1; icol <= ncol; icol++) {
                    var indx = icol-1 + (irow-1)*pmtr[ipmtr].ncol;

                    newTD = document.createElement("TD");
                    newTR.appendChild(newTD);

                    var newInput = document.createElement("input");
                    newInput.type  = "text";
                    newInput.name  = "row"+irow+"col"+icol+"val";
                    newInput.size  = 12;
                    newInput.value = pmtr[ipmtr].value[indx];
                    newTD.appendChild(newInput);

                    if (irow == 1 && icol == 1) {
                        wv.getFocus = newInput;
                    }
                }
            }
        }

        if (pmtr[ipmtr].type == OCSM_DESPMTR) {
            document.getElementById("editValuDot").hidden = false;

            var valuDotTable = document.getElementById("editValuDotTable");

            // remove old table entries
            if (valuDotTable) {
                var child1 = valuDotTable.lastChild;
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
                    valuDotTable.removeChild(child1);
                    child1 = valuDotTable.lastChild;
                }
            }

            // build the table that will contain velocities
            for (var irow = 0; irow <= nrow; irow++) {
                var newTR = document.createElement("TR");
                valuDotTable.appendChild(newTR);

                // fill the row
                if (irow == 0) {
                    var newTD = document.createElement("TD");
                    newTR.appendChild(newTD);

                    var newText = document.createTextNode("");
                    newTD.appendChild(newText);

                    for (var icol = 1; icol <= ncol; icol++) {
                        newTD = document.createElement("TD");
                        newTR.appendChild(newTD);

                        newText = document.createTextNode("col\u00a0"+icol);
                        newTD.appendChild(newText);
                    }
                } else{
                    var newTD = document.createElement("TD");
                    newTR.appendChild(newTD);

                    var newText = document.createTextNode("row\u00a0"+irow);
                    newTD.appendChild(newText);

                    for (var icol = 1; icol <= ncol; icol++) {
                        var indx = icol-1 + (irow-1)*pmtr[ipmtr].ncol;

                        newTD = document.createElement("TD");
                        newTR.appendChild(newTD);

                        var newInput = document.createElement("input");
                        newInput.type  = "text";
                        newInput.name  = "row"+irow+"col"+icol+"dot";
                        newInput.size  = 12;
                        newInput.value = pmtr[ipmtr].dot[indx];
                        newTD.appendChild(newInput);
                    }
                }
            }
        } else {
            document.getElementById("editValuDot").hidden = true;
        }

    // setting up for a Caps Value
    } else if (wv.curCval >= 0) {
        var icval = wv.curCval;
        var name  = cval[icval].name;
        var nrow  = cval[icval].nrow;
        var ncol  = cval[icval].ncol;

        var editValuForm = document.getElementById("editValuForm");

        // fill in the Parameter name
        document.getElementById("valuName").firstChild["data"] = name;

        var valuValTable = document.getElementById("editValuValTable");

        // remove old table entries
        if (valuValTable) {
            var child1 = valuValTable.lastChild;
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
                valuValTable.removeChild(child1);
                child1 = valuValTable.lastChild;
            }
        }

        // build the table that will contain values
        for (var irow = 0; irow <= nrow; irow++) {
            var newTR = document.createElement("TR");
            valuValTable.appendChild(newTR);

            // fill the row
            if (irow == 0) {
                var newTD = document.createElement("TD");
                newTR.appendChild(newTD);

                var newText = document.createTextNode("");
                newTD.appendChild(newText);

                for (var icol = 1; icol <= ncol; icol++) {
                    newTD = document.createElement("TD");
                    newTR.appendChild(newTD);

                    newText = document.createTextNode("col\u00a0"+icol);
                    newTD.appendChild(newText);
                }
            } else{
                var newTD = document.createElement("TD");
                newTR.appendChild(newTD);

                var newText = document.createTextNode("row\u00a0"+irow);
                newTD.appendChild(newText);

                for (var icol = 1; icol <= ncol; icol++) {
                    var indx = icol-1 + (irow-1)*cval[icval].ncol;

                    newTD = document.createElement("TD");
                    newTR.appendChild(newTD);

                    var newInput = document.createElement("input");
                    newInput.type  = "text";
                    newInput.name  = "row"+irow+"col"+icol+"val";
                    newInput.size  = 12;
                    newInput.value = cval[icval].value[indx];
                    newTD.appendChild(newInput);

                    if (irow == 1 && icol == 1) {
                        wv.getFocus = newInput;
                    }
                }
            }
        }

        document.getElementById("editValuDot").hidden = true;

    // ooops, problem
    } else {
        alert("we have a problem: curPmtr="+wv.curPmtr+" curCval="+wv.curCval);
    }

    return 0;
};


//
// send a message to the server
//
var browserToServer = function (text) {
    // alert("in browserToServer(text="+text+")");

    if (wv.debugUI) {
        var date = new Date;
        console.log("("+date.toTimeString().substring(0,8)+") browser-->server: "+text.substring(0,40));
    }

    wv.socketUt.send(text);
};


//
// count number of Parameter value changes
//
var numberOfPmtrChanges = function () {
    // alert("in numberOfPmtrChanges()");

    var nchange = 0;

    var editValuForm = document.getElementById("editValuForm");

    var ipmtr = wv.curPmtr;
    var irow;
    var icol;

    // examine each of the values
    var index   = -1;
    for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
        for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
            index++;

            // get the new value
            var myInput = editValuForm["row"+irow+"col"+icol+"val"];
            var value   = myInput.value.replace(/\s/g, "");

            if (value != pmtr[ipmtr].value[index]) {
                nchange++;
            }
        }
    }

    // examine each of the velocities
    if (pmtr[ipmtr].type == OCSM_DESPMTR) {
        index   = -1;
        for (irow = 1; irow <= pmtr[ipmtr].nrow; irow++) {
            for (icol = 1; icol <= pmtr[ipmtr].ncol; icol++) {
                index++;

                // get the new value
                var myInput = editValuForm["row"+irow+"col"+icol+"dot"];
                var value   = myInput.value.replace(/\s/g, "");

                if (value != pmtr[ipmtr].dot[index]) {
                    nchange++;
                }
            }
        }
    }

    // return the number of changes
    return nchange;
};


//
// count number of Branch changes
//
var numberOfBrchChanges = function () {
    // alert("in numberOfBrchChanges()");

    var nchange = 0;

    var editBrchForm = document.getElementById("editBrchForm");

    var ibrch = wv.curBrch;

    // check the name
    if (brch[ibrch].name != editBrchForm.brchName.value.replace(/\s/g, "")) {
        nchange++;
    }

    // check the activity
    if (brch[ibrch].actv != editBrchForm.activity.value) {
        nchange++;
    }

    // check the arguments
    for (var iarg = 0; iarg < wv.numArgs; iarg++) {
        if        (iarg == 0) {
            var name  = document.getElementById("argName1").firstChild["data"];
            var value = editBrchForm.            argValu1.value.replace(/\s/g, "");
        } else if (iarg == 1) {
            var name  = document.getElementById("argName2").firstChild["data"];
            var value = editBrchForm.            argValu2.value.replace(/\s/g, "");
        } else if (iarg == 2) {
            var name  = document.getElementById("argName3").firstChild["data"];
            var value = editBrchForm.            argValu3.value.replace(/\s/g, "");
        } else if (iarg == 3) {
            var name  = document.getElementById("argName4").firstChild["data"];
            var value = editBrchForm.            argValu4.value.replace(/\s/g, "");
        } else if (iarg == 4) {
            var name  = document.getElementById("argName5").firstChild["data"];
            var value = editBrchForm.            argValu5.value.replace(/\s/g, "");
        } else if (iarg == 5) {
            var name  = document.getElementById("argName6").firstChild["data"];
            var value = editBrchForm.            argValu6.value.replace(/\s/g, "");
        } else if (iarg == 6) {
            var name  = document.getElementById("argName7").firstChild["data"];
            var value = editBrchForm.            argValu7.value.replace(/\s/g, "");
        } else if (iarg == 7) {
            var name  = document.getElementById("argName8").firstChild["data"];
            var value = editBrchForm.            argValu8.value.replace(/\s/g, "");
        } else if (iarg == 8) {
            var name  = document.getElementById("argName9").firstChild["data"];
            var value = editBrchForm.            argValu9.value.replace(/\s/g, "");
        }

        var output;
        if (name.charAt(0) != "$" || value.length <= 0) {
            output =       value;
        } else {
            output = "$" + value;
        }

        if (output != brch[ibrch].args[iarg]) {
            nchange++;
        }
    }

    // check the attributes
    for (var iattr = 0; iattr < brch[ibrch].attrs.length; iattr++) {
        if        (iattr == 0) {
            var name  = document.getElementById("attrName1").firstChild["data"];
            var value = editBrchForm.            attrValu1.value.replace(/\s/g, "");
        } else if (iattr == 1) {
            var name  = document.getElementById("attrName2").firstChild["data"];
            var value = editBrchForm.            attrValu2.value.replace(/\s/g, "");
        } else if (iattr == 2) {
            var name  = document.getElementById("attrName3").firstChild["data"];
            var value = editBrchForm.            attrValu3.value.replace(/\s/g, "");
        } else if (iattr == 3) {
            var name  = document.getElementById("attrName4").firstChild["data"];
            var value = editBrchForm.            attrValu4.value.replace(/\s/g, "");
        } else if (iattr == 4) {
            var name  = document.getElementById("attrName5").firstChild["data"];
            var value = editBrchForm.            attrValu5.value.replace(/\s/g, "");
        } else if (iattr == 5) {
            var name  = document.getElementById("attrName6").firstChild["data"];
            var value = editBrchForm.            attrValu6.value.replace(/\s/g, "");
        } else if (iattr == 6) {
            var name  = document.getElementById("attrName7").firstChild["data"];
            var value = editBrchForm.            attrValu7.value.replace(/\s/g, "");
        } else if (iattr == 7) {
            var name  = document.getElementById("attrName8").firstChild["data"];
            var value = editBrchForm.            attrValu8.value.replace(/\s/g, "");
        } else if (iattr == 8) {
            var name  = document.getElementById("attrName9").firstChild["data"];
            var value = editBrchForm.            attrValu9.value.replace(/\s/g, "");
        }

        if (value != brch[ibrch].attrs[iattr][2]) {
            nchange++;
        }
    }

    // return the number of changes
    return nchange;
};


//
// unhighlight background colors in column 1 of myTree
//
var unhighlightColumn1 = function () {
    // alert("in unhighlightColumn1()");

    for (var jnode = 0; jnode < myTree.name.length; jnode++) {
        var myElem = document.getElementById("node"+jnode+"col1");
        if (myElem) {
            if (myElem.className == "currentTD" || myElem.className == "parentTD" || myElem.className == "childTD") {
                myElem.className = undefined;
            }
        }
    }
};


//
// change a property in the TreeWindow
//
var changeProp = function (inode, icol, state) {
    // alert("in changeProp(inode="+inode+", icol="+icol+", state="+state+")");

    myTree.prop(inode, icol, state);

    if (wv.myRole == 0 && wv.numUsers > 1) {
        browserToServer("toggle|"+inode+"|"+icol+"|"+state+"|");
    }
};


//
// callback from "Copy" button in .csm editor
//
var cmdEditCopy = function (cm) {
    // alert("in cmdEditCopy(cm="+cm"+)");

    // copy to private clipboard
    wv.clipboard = cm.getRange(cm.getCursor("start"), cm.getCursor("end"));

    // focus back to textarea
    cm.focus();
};


//
// callback from "Cut" button in .csm editor
//
var cmdEditCut = function (cm) {
    // alert("in cmdEditCut(cm="+cm+")");

    if (checkIfWithBall() === false) {
        return;
    }

    // copy to private clipboard
    wv.clipboard = cm.getRange(cm.getCursor("start"), cm.getCursor("end"));

    // remove from textarea
    cm.replaceSelection("", null, "paste");

    // focus back to textarea
    cm.focus();
};


//
// callback from "Paste" button in .csm editor
//
var cmdEditPaste = function (cm) {
    // alert("in cmdEditPaste(cm="+cm+")");

    if (checkIfWithBall() === false) {
        return;
    }

    // copy from private clipboard to textarea
    cm.replaceSelection(wv.clipboard, null, "paste")

    // focus back to textarea
    cm.focus();
};


//
// callback from "Search" button in .csm editor
var cmdEditFind = function (cm) {
    // alert("in cmdEditFind(cm="+cm+")");

    CodeMirror.commands.find(cm);
};


//
// callback from "Next" button in .csm editor
//
var cmdEditFindNext = function (cm) {
    // alert("in cmdEditFindNext(cm="+cm+")");

    CodeMirror.commands.findNext(cm);
};


//
// callback from "Prev" button in .csm editor
//
var cmdEditFindPrev = function (cm) {
    // alert("in cmdEditFindPrev(cm="+cm+")");

    CodeMirror.commands.findPrev(cm);
};


//
// callback from "Replace" button in .csm editor
//
var cmdEditReplace = function (cm) {
    // alert("in cmdEditReplace(cm="+cm+")");

    if (checkIfWithBall() === false) {
        return;
    }

    CodeMirror.commands.replace(cm);
};


//
// callback from "Comment" button in .csm editor
//
var cmdEditComment = function () {
    // alert("in cmdEditComment()");

    if (checkIfWithBall() === false) {
        return;
    }

    // find lines in region (or line with current cursor)
    var begLine = wv.codeMirror.getCursor("start").line;
    var endLine = wv.codeMirror.getCursor("end").line;

    // do not include line at end if cursor is at beginning of line
    if (wv.codeMirror.getCursor("end").ch == 0 && endLine > begLine) {
        endLine--;
    }

    // determine if we are adding or removing comment based upon first line
    if (wv.codeMirror.getLine(begLine).match(/^#--- /) !== null) {

        // remove comments
        for (var iline = begLine; iline <= endLine; iline++) {
            if (wv.codeMirror.getLine(iline).match(/^#--- /) !== null) {
                wv.codeMirror.replaceRange("", CodeMirror.Pos(iline,0), CodeMirror.Pos(iline,5));
            }
        }
    } else {

        // add comments
        for (var iline = begLine; iline <= endLine; iline++) {
            wv.codeMirror.replaceRange("#--- ", CodeMirror.Pos(iline,0), CodeMirror.Pos(iline,0));
        }
    }
};


//
// callback from "Indent" button in .csm editor
//
var cmdEditIndent = function () {
    // alert("in cmdIndent()");

    if (checkIfWithBall() === false) {
        return;
    }

    // find lines in region (or line with current cursor)
    var begLine = wv.codeMirror.getCursor("start").line;
    var endLine = wv.codeMirror.getCursor("end").line;

    // do not include line at end if cursor is at beginning of line
    if (wv.codeMirror.getCursor("end").ch == 0 && endLine > begLine) {
        endLine--;
    }

    // auto indent
    for (var iline = begLine; iline <= endLine; iline++) {
        wv.codeMirror.indentLine(iline, "smart");

        // special indentation for "attribute", "csystem", and "name"
        var ichar = wv.codeMirror.getLine(iline).match(/^\s*(attribute|ATTRIBUTE|csystem|CSYSTEM|name|NAME)\s/);

        if (ichar !== null) {
            wv.codeMirror.replaceRange("   ", CodeMirror.Pos(iline,ichar), CodeMirror.Pos(iline,ichar));
        }
    }
};


//
// callback from "Hint" button in .csm editor
//
var cmdEditHint = function () {
    // alert("in cmdEditHint");

    // get the current line
    var curLine = wv.codeMirror.getLine(wv.codeMirror.getCursor("start").line);

    // check the current line type
    var hintText = "*** no hint available for: " + curLine;

    if        (curLine.match(/^\s*$/) !== null) {
        hintText = "*** no hint available for blank line";
    } else if (curLine.match(/^\s*#/) !== null) {
        hintText =        "hint:: # comment";
    } else if (curLine.match(/^\s*applycsys/i) !== null) {
        hintText =        "hint:: APPLYCSYS $csysName ibody=0";
    } else if (curLine.match(/^\s*arc/i) !== null) {
        hintText =        "hint:: ARC       xend yend zend dist $plane=xy";
    } else if (curLine.match(/^\s*assert/i) !== null) {
        hintText =        "hint:: ASSERT    arg1 arg2 toler=0 verify=0";
    } else if (curLine.match(/^\s*attribute/i) !== null) {
        hintText =        "hint:: ATTRIBUTE $attrName attrValue";
    } else if (curLine.match(/^\s*bezier/i) !== null) {
        hintText =        "hint:: BEZIER    x y z";
    } else if (curLine.match(/^\s*blend/i) !== null) {
        hintText =        "hint:: BLEND     begList=0 endList=0 reorder=0 oneFace=0 periodic=0";
    } else if (curLine.match(/^\s*box/i) !== null) {
        hintText =        "hint:: BOX       xbase ybase zbase dx dy dz";
    } else if (curLine.match(/^\s*catbeg/i) !== null) {
        hintText =        "hint:: CATBEG    sigCode";
    } else if (curLine.match(/^\s*catend/i) !== null) {
        hintText =        "hint:: CATEND";
    } else if (curLine.match(/^\s*cfgpmtr/i) !== null) {
        hintText =        "hint:: CFGPMTR   $pmtrName values";
    } else if (curLine.match(/^\s*chamfer/i) !== null) {
        hintText =        "hint:: CHAMFER   radius edgeList=0 listStyle=0";
    } else if (curLine.match(/^\s*cirarc/i) !== null) {
        hintText =        "hint:: CIRARC    xon yon zon xend yend zend";
    } else if (curLine.match(/^\s*combine/i) !== null) {
        hintText =        "hint:: COMBINE   toler=0";
    } else if (curLine.match(/^\s*cone/i) !== null) {
        hintText =        "hint:: CONE      xvrtx yvrtx zvrtx xbase ybase zbase radius";
    } else if (curLine.match(/^\s*connect/i) !== null) {
        hintText =        "hint:: CONNECT   faceList1 faceList2 edgeList1=0 edgeList2=0 toler=0";
    } else if (curLine.match(/^\s*conpmtr/i) !== null) {
        hintText =        "hint:: CONPMTR   $pmtrName value";
    } else if (curLine.match(/^\s*csystem/i) !== null) {
        hintText =        "hint:: CSYSTEM   $csysName csysList";
    } else if (curLine.match(/^\s*cylinder/i) !== null) {
        hintText =        "hint:: CYLINDER  xbeg ybeg zbeg xend yend zend radius";
    } else if (curLine.match(/^\s*despmtr/i) !== null) {
        hintText =        "hint:: DESPMTR   $pmtrName values";
    } else if (curLine.match(/^\s*dimension/i) !== null) {
        hintText =        "hint:: DIMENSION $pmtrName nrow ncol despmtr=0";
    } else if (curLine.match(/^\s*dump/i) !== null) {
        hintText =        "hint:: DUMP      $filename remove=0 toMark=0 withTess=0";
    } else if (curLine.match(/^\s*elseif/i) !== null) {
        hintText =        "hint:: ELSEIF    val1 $op1 val2 $op2=and val3=0 $op3=eq val4=0";
    } else if (curLine.match(/^\s*else/i) !== null) {
        hintText =        "hint:: ELSE";
    } else if (curLine.match(/^\s*endif/i) !== null) {
        hintText =        "hint:: ENDIF";
    } else if (curLine.match(/^\s*end/i) !== null) {
        hintText =        "hint:: END";
    } else if (curLine.match(/^\s*evaluate/i) !== null) {
        hintText =        "hint:: EVALUATE  $type arg1 ...";
    } else if (curLine.match(/^\s*extract/i) !== null) {
        hintText =        "hint:: EXTRACT   index";
    } else if (curLine.match(/^\s*extrude/i) !== null) {
        hintText =        "hint:: EXTRUDE   dx dy dz";
    } else if (curLine.match(/^\s*fillet/i) !== null) {
        hintText =        "hint:: FILLET    radius edgeList=0 listStyle=0";
    } else if (curLine.match(/^\s*getattr/i) !== null) {
        hintText =        "hint:: GETATTR   $pmtrName attrID";
    } else if (curLine.match(/^\s*group/i) !== null) {
        hintText =        "hint:: GROUP     nbody=0";
    } else if (curLine.match(/^\s*hollow/i) !== null) {
        hintText =        "hint:: HOLLOW    thick=0 entList=0 listStyle=0";
    } else if (curLine.match(/^\s*ifthen/i) !== null) {
        hintText =        "hint:: IFTHEN    val1 $op1 val2 $op2=and val3=0 $op3=eq val4=0";
    } else if (curLine.match(/^\s*import/i) !== null) {
        hintText =        "hint:: IMPORT    $filename bodynumber=1";
    } else if (curLine.match(/^\s*interface/i) !== null) {
        hintText =        "hint:: INTERFACE $argName $argType default";
    } else if (curLine.match(/^\s*intersect/i) !== null) {
        hintText =        "hint:: INTERSECT $order=none index=1 maxtol=0";
    } else if (curLine.match(/^\s*join/i) !== null) {
        hintText =        "hint:: JOIN      toler=0";
    } else if (curLine.match(/^\s*lbound/i) !== null) {
        hintText =        "hint:: LBOUND    $pmtrName bounds";
    } else if (curLine.match(/^\s*linseg/i) !== null) {
        hintText =        "hint:: LINSEG    x y z";
    } else if (curLine.match(/^\s*loft/i) !== null) {
        hintText =        "hint:: LOFT      smooth (DEPRECATED)";
    } else if (curLine.match(/^\s*macbeg/i) !== null) {
        hintText =        "hint:: MACBEG    imacro (DEPRECATED)";
    } else if (curLine.match(/^\s*macend/i) !== null) {
        hintText =        "hint:: MACEND (deprecate)";
    } else if (curLine.match(/^\s*mark/i) !== null) {
        hintText =        "hint:: MARK";
    } else if (curLine.match(/^\s*message/i) !== null) {
        hintText =        "hint:: MESSAGE   $text $schar=_";
    } else if (curLine.match(/^\s*mirror/i) !== null) {
        hintText =        "hint:: MIRROR    nx ny nz dist=0";
    } else if (curLine.match(/^\s*name/i) !== null) {
        hintText =        "hint:: NAME      $branchName";
    } else if (curLine.match(/^\s*outpmtr/i) !== null) {
        hintText =        "hint:: OUTPMTR   $pmtrName";
    } else if (curLine.match(/^\s*patbeg/i) !== null) {
        hintText =        "hint:: PATBEG    $pmtrName ncopy";
    } else if (curLine.match(/^\s*patbreak/i) !== null) {
        hintText =        "hint:: PATBREAK  expr";
    } else if (curLine.match(/^\s*patend/i) !== null) {
        hintText =        "hint:: PATEND";
    } else if (curLine.match(/^\s*point/i) !== null) {
        hintText =        "hint:: POINT     xloc yloc zloc";
    } else if (curLine.match(/^\s*project/i) !== null) {
        hintText =        "hint:: PROJECT   x y z dx dy dz useEdges=0";
    } else if (curLine.match(/^\s*recall/i) !== null) {
        hintText =        "hint:: RECALL    imacro (DEPRECATED)";
    } else if (curLine.match(/^\s*reorder/i) !== null) {
        hintText =        "hint:: REORDER   ishift iflip=0";
    } else if (curLine.match(/^\s*restore/i) !== null) {
        hintText =        "hint:: RESTORE   $name index=0";
    } else if (curLine.match(/^\s*revolve/i) !== null) {
        hintText =        "hint:: REVOLVE   xorig yorig zorig dxaxis dyaxis dzaxis angDeg";
    } else if (curLine.match(/^\s*rotatex/i) !== null) {
        hintText =        "hint:: ROTATEX   angDeg yaxis=0 zaxis=0";
    } else if (curLine.match(/^\s*rotatey/i) !== null) {
        hintText =        "hint:: ROTATEY   angDeg zaxis=0 xaxis=0";
    } else if (curLine.match(/^\s*rotatez/i) !== null) {
        hintText =        "hint:: ROTATEZ   angDeg xaxis=0 yaxis=0";
    } else if (curLine.match(/^\s*rule/i) !== null) {
        hintText =        "hint:: RULE      reorder=0 periodic=0";
    } else if (curLine.match(/^\s*scale/i) !== null) {
        hintText =        "hint:: SCALE     fact xcent=0 ycent=0 zcent=0";
    } else if (curLine.match(/^\s*select/i) !== null) {
        hintText =        "hint:: SELECT    $type arg1 ...";
    } else if (curLine.match(/^\s*set/i) !== null) {
        hintText =        "hint:: SET       $pmtrName exprs";
    } else if (curLine.match(/^\s*skbeg/i) !== null) {
        hintText =        "hint:: SKBEG     x y z relative=0";
    } else if (curLine.match(/^\s*skcon/i) !== null) {
        hintText =        "hint:: SKCON     $type index1 index2=-1 $value=0";
    } else if (curLine.match(/^\s*skend/i) !== null) {
        hintText =        "hint:: SKEND     wireonly=0";
    } else if (curLine.match(/^\s*skvar/i) !== null) {
        hintText =        "hint:: SKVAR     $type valList";
    } else if (curLine.match(/^\s*solbeg/i) !== null) {
        hintText =        "hint:: SOLBEG    $varList";
    } else if (curLine.match(/^\s*solcon/i) !== null) {
        hintText =        "hint:: SOLCON    $expr";
    } else if (curLine.match(/^\s*solend/i) !== null) {
        hintText =        "hint:: SOLEND";
    } else if (curLine.match(/^\s*sphere/i) !== null) {
        hintText =        "hint:: SPHERE    xcent ycent zcent radius";
    } else if (curLine.match(/^\s*spline/i) !== null) {
        hintText =        "hint:: SPLINE    x y z";
    } else if (curLine.match(/^\s*sslope/i) !== null) {
        hintText =        "hint:: SSLOPE    dx dy dz";
    } else if (curLine.match(/^\s*store/i) !== null) {
        hintText =        "hint:: STORE     $name index=0 keep=0";
    } else if (curLine.match(/^\s*subtract/i) !== null) {
        hintText =        "hint:: SUBTRACT  $order=none index=1 maxtol=0";
    } else if (curLine.match(/^\s*sweep/i) !== null) {
        hintText =        "hint:: SWEEP";
    } else if (curLine.match(/^\s*throw/i) !== null) {
        hintText =        "hint:: THROW     sigCode";
    } else if (curLine.match(/^\s*torus/i) !== null) {
        hintText =        "hint:: TORUS     xcent ycent zcent dxaxis dyaxis dzaxis majorRad minorRad";
    } else if (curLine.match(/^\s*translate/i) !== null) {
        hintText =        "hint:: TRANSLATE dx dy dz";
    } else if (curLine.match(/^\s*ubound/i) !== null) {
        hintText =        "hint:: UBOUND    $pmtrName bounds";
    } else if (curLine.match(/^\s*udparg/i) !== null) {
        hintText =        "hint:: UDPARG    $primtype $argName1 argValue1 $argName2 argValue2 $argName3 argValue3 $argName4 argValue4";
    } else if (curLine.match(/^\s*udprim/i) !== null) {
        hintText =        "hint:: UDPRIM    $primtype $argName1 argValue1 $argName2 argValue2 $argName3 argValue3 $argName4 argValue4";
    } else if (curLine.match(/^\s*union/i) !== null) {
        hintText =        "hint:: UNION     toMark=0 trimList=0 maxtol=0";
    }

    // post the hint (for at least 15 seconds)
    wv.codeMirror.openNotification(hintText, {duration: 15000});

    // focus back to textarea
    wv.codeMirror.focus();
};


//
// callback from "Undo" button in .csm editor
//
var cmdEditUndo = function (cm) {
    // alert("in cmdEditUndo()");

    if (checkIfWithBall() === false) {
        return;
    }

    CodeMirror.commands.undo(cm);
};


//
// callback from "Debug" button in .csm editor
//
var cmdEditDebug = function (cm) {
    // alert("in cmdEditDebug()");

    if (checkIfWithBall() === false) {
        return;
    }

    // get the current line
    var curLine = wv.codeMirror.getCursor("start").line + 1;

    // get filename
    var filelist = wv.filenames.split("|");
    if (wv.fileindx === undefined && filelist.length == 3) {
        wv.fileindx = 1;
    }

    // get the info from the server
    browserToServer("getBodyDetails|"+filelist[wv.fileindx]+"|"+curLine+"|");
};


//
// print an object and its contents
//
var printObject = function (obj) {
    var out = '';

    for (var p in obj) {
        out += p + ': ' + obj[p] + '\n';
    }

    alert(out);
};


//
// simple text formatter patterned after C's sprintf
//
var sprintf = function () {

    // if no arguments, return an empty string
    var narg = arguments.length;
    if (narg == 0) {
        return "";
    }

    // otherwise, build output from input
    var format = arguments[0];
    var answer = "";
    var iarg   = 1;

    while (1) {
        var indx = format.indexOf("%");
        if (indx >= 0) {
            answer += format.substring(0, indx);
            if (iarg < narg) {
                answer += arguments[iarg++];
            } else {
                answer += "%";
            }
            format  = format.substring(indx+1);
        } else {
            answer += format;
            break;
        }
    }

    return answer;
};


//
// definition of the ".csm" mode to be used by CodeMirror
//
CodeMirror.defineSimpleMode("csm_mode", {

  // The start state contains the rules that are intially used
  start: [

    {token: "comment", regex: /#.*/},
    {token: "comment", regex: /\\.*/},

    {token: "string",  regex: /\$[^+,')# ]*('[+,)][^+,)'# ]*)*/},

    {token: "keyword", regex: /\s*(applycsys|APPLYCSYS)\b/, sol: true},
    {token: "keyword", regex: /\s*(arc|ARC)\b/,             sol: true},
    {token: "keyword", regex: /\s*(assert|ASSERT)\b/,       sol: true},
    {token: "keyword", regex: /\s*(attribute|ATTRIBUTE)\b/, sol: true},
    {token: "keyword", regex: /\s*(bezier|BEZIER)\b/,       sol: true},
    {token: "keyword", regex: /\s*(blend|BLEND)\b/,         sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(box|BOX)\b/,             sol: true},
    {token: "keyword", regex: /\s*(catbeg|CATBEG)\b/,       sol: true, indent: true},
    {token: "keyword", regex: /\s*(catend|CATEND)\b/,       sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(cfgpmtr|CFGPMTR)\b/,     sol: true},
    {token: "keyword", regex: /\s*(chamfer|CHAMFER)\b/,     sol: true},
    {token: "keyword", regex: /\s*(cirarc|CIRARC)\b/,       sol: true},
    {token: "keyword", regex: /\s*(combine|COMBINE)\b/,     sol: true},
    {token: "keyword", regex: /\s*(cone|CONE)\b/,           sol: true},
    {token: "keyword", regex: /\s*(connect|CONNECT)\b/,     sol: true},
    {token: "keyword", regex: /\s*(conpmtr|CONPMTR)\b/,     sol: true},
    {token: "keyword", regex: /\s*(csystem|CSYSTEM)\b/,     sol: true},
    {token: "keyword", regex: /\s*(cylinder|CYLINDER)\b/,   sol: true},
    {token: "keyword", regex: /\s*(despmtr|DESPMTR)\b/,     sol: true},
    {token: "keyword", regex: /\s*(dimension|DIMENSION)\b/, sol: true},
    {token: "keyword", regex: /\s*(dump|DUMP)\b/,           sol: true},
    {token: "keyword", regex: /\s*(else|ELSE)\b/,           sol: true, indent: true, dedent: true},
    {token: "keyword", regex: /\s*(elseif|ELSEIF)\b/,       sol: true, indent: true, dedent: true},
    {token: "keyword", regex: /\s*(end|END)\b/,             sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(endif|ENDIF)\b/,         sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(evaluate|EVALUATE)\b/,   sol: true},
    {token: "keyword", regex: /\s*(extract|EXTRACT)\b/,     sol: true},
    {token: "keyword", regex: /\s*(extrude|EXTRUDE)\b/,     sol: true},
    {token: "keyword", regex: /\s*(fillet|FILLET)\b/,       sol: true},
    {token: "keyword", regex: /\s*(getattr|GETATTR)\b/,     sol: true},
    {token: "keyword", regex: /\s*(group|GROUP)\b/,         sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(hollow|HOLLOW)\b/,       sol: true},
    {token: "keyword", regex: /\s*(ifthen|IFTHEN)\b/,       sol: true, indent: true},
    {token: "keyword", regex: /\s*(import|IMPORT)\b/,       sol: true},
    {token: "keyword", regex: /\s*(interface|INTERFACE)\b/, sol: true},
    {token: "keyword", regex: /\s*(intersect|INTERSECT)\b/, sol: true},
    {token: "keyword", regex: /\s*(join|JOIN)\b/,           sol: true},
    {token: "keyword", regex: /\s*(lbound|LBOUND)\b/,       sol: true},
    {token: "keyword", regex: /\s*(linseg|LINSEG)\b/,       sol: true},
    {token: "keyword", regex: /\s*(loft|LOFT)\b/,           sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(macbeg|MACBEG)\b/,       sol: true, indent: true},
    {token: "keyword", regex: /\s*(macend|MACEND)\b/,       sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(mark|MARK)\b/,           sol: true, indent: true},
    {token: "keyword", regex: /\s*(message|MESSAGE)\b/,     sol: true},
    {token: "keyword", regex: /\s*(mirror|MIRROR)\b/,       sol: true},
    {token: "keyword", regex: /\s*(name|NAME)\b/,           sol: true},
    {token: "keyword", regex: /\s*(outpmtr|OUTPMTR)\b/,     sol: true},
    {token: "keyword", regex: /\s*(patbeg|PATBEG)\b/,       sol: true, indent: true},
    {token: "keyword", regex: /\s*(patbreak|PATBREAK)\b/,   sol: true, indent: true, dedent: true},
    {token: "keyword", regex: /\s*(patend|PATEND)\b/,       sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(point|POINT)\b/,         sol: true},
    {token: "keyword", regex: /\s*(project|PROJECT)\b/,     sol: true},
    {token: "keyword", regex: /\s*(recall|RECALL)\b/,       sol: true},
    {token: "keyword", regex: /\s*(reorder|REORDER)\b/,     sol: true},
    {token: "keyword", regex: /\s*(restore|RESTORE)\b/,     sol: true},
    {token: "keyword", regex: /\s*(revolve|REVOLVE)\b/,     sol: true},
    {token: "keyword", regex: /\s*(rotatex|ROTATEX)\b/,     sol: true},
    {token: "keyword", regex: /\s*(rotatey|ROTATEY)\b/,     sol: true},
    {token: "keyword", regex: /\s*(rotatez|ROTATEZ)\b/,     sol: true},
    {token: "keyword", regex: /\s*(rule|RULE)\b/,           sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(scale|SCALE)\b/,         sol: true},
    {token: "keyword", regex: /\s*(select|SELECT)\b/,       sol: true},
    {token: "keyword", regex: /\s*(set|SET)\b/,             sol: true},
    {token: "keyword", regex: /\s*(skbeg|SKBEG)\b/,         sol: true, indent: true},
    {token: "keyword", regex: /\s*(skcon|SKCON)\b/,         sol: true},
    {token: "keyword", regex: /\s*(skend|SKEND)\b/,         sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(skvar|SKVAR)\b/,         sol: true},
    {token: "keyword", regex: /\s*(solbeg|SOLBEG)\b/,       sol: true, indent: true},
    {token: "keyword", regex: /\s*(solcon|SOLCON)\b/,       sol: true},
    {token: "keyword", regex: /\s*(solend|SOLEND)\b/,       sol: true,               dedent: true},
    {token: "keyword", regex: /\s*(sphere|SPHERE)\b/,       sol: true},
    {token: "keyword", regex: /\s*(spline|SPLINE)\b/,       sol: true},
    {token: "keyword", regex: /\s*(sslope|SSLOPE)\b/,       sol: true},
    {token: "keyword", regex: /\s*(store|STORE)\b/,         sol: true},
    {token: "keyword", regex: /\s*(subtract|SUBTRACT)\b/,   sol: true},
    {token: "keyword", regex: /\s*(sweep|SWEEP)\b/,         sol: true},
    {token: "keyword", regex: /\s*(throw|THROW)\b/,         sol: true},
    {token: "keyword", regex: /\s*(torus|TORUS)\b/,         sol: true},
    {token: "keyword", regex: /\s*(translate|TRANSLATE)\b/, sol: true},
    {token: "keyword", regex: /\s*(ubound|UBOUND)\b/,       sol: true},
    {token: "keyword", regex: /\s*(udparg|UDPARG)\b/,       sol: true},
    {token: "keyword", regex: /\s*(udprim|UDPRIM)\b/,       sol: true},
    {token: "keyword", regex: /\s*(union|UNION)\b/,         sol: true},

    // positive lookaround to make sure followed by (
    {token: "builtin", regex: /pi(?=\()/},
    {token: "builtin", regex: /min(?=\()/},
    {token: "builtin", regex: /max(?=\()/},
    {token: "builtin", regex: /sqrt(?=\()/},
    {token: "builtin", regex: /abs(?=\()/},
    {token: "builtin", regex: /int(?=\()/},
    {token: "builtin", regex: /nint(?=\()/},
    {token: "builtin", regex: /ceil(?=\()/},
    {token: "builtin", regex: /floor(?=\()/},
    {token: "builtin", regex: /mod(?=\()/},
    {token: "builtin", regex: /sign(?=\()/},
    {token: "builtin", regex: /exp(?=\()/},
    {token: "builtin", regex: /log(?=\()/},
    {token: "builtin", regex: /log10(?=\()/},
    {token: "builtin", regex: /sind?(?=\()/},
    {token: "builtin", regex: /asind?(?=\()/},
    {token: "builtin", regex: /cosd?(?=\()/},
    {token: "builtin", regex: /acosd?(?=\()/},
    {token: "builtin", regex: /tand?(?=\()/},
    {token: "builtin", regex: /atan2?d?(?=\()/},
    {token: "builtin", regex: /hypot3?(?=\()/},
    {token: "builtin", regex: /incline(?=\()/},
    {token: "builtin", regex: /[XY]cent(?=\()/},
    {token: "builtin", regex: /[XY]midl(?=\()/},
    {token: "builtin", regex: /seglen(?=\()/},
    {token: "builtin", regex: /radius(?=\()/},
    {token: "builtin", regex: /sweep(?=\()/},
    {token: "builtin", regex: /turnang(?=\()/},
    {token: "builtin", regex: /dip(?=\()/},
    {token: "builtin", regex: /smallang(?=\()/},
    {token: "builtin", regex: /val2str(?=\()/},
    {token: "builtin", regex: /str2val(?=\()/},
    {token: "builtin", regex: /findstr(?=\()/},
    {token: "builtin", regex: /slice(?=\()/},
    {token: "builtin", regex: /path(?=\()/},
    {token: "builtin", regex: /ifzero(?=\()/},
    {token: "builtin", regex: /ifpos(?=\()/},
    {token: "builtin", regex: /ifneg(?=\()/},
    {token: "builtin", regex: /ifmatch(?=\()/},
    {token: "builtin", regex: /ifnan(?=\()/},

    // negative lookaround to make sure not followed by letter
    {token: "meta", regex: /\.nrow(?!\w)/},
    {token: "meta", regex: /\.ncol(?!\w)/},
    {token: "meta", regex: /\.size(?!\w)/},
    {token: "meta", regex: /\.sum(?!\w)/},
    {token: "meta", regex: /\.norm(?!\w)/},
    {token: "meta", regex: /\.min(?!\w)/},
    {token: "meta", regex: /\.max(?!\w)/},

    {token: "special", regex: /\s(eq|EQ)\s/},
    {token: "special", regex: /\s(ne|NE)\s/},
    {token: "special", regex: /\s(lt|LT)\s/},
    {token: "special", regex: /\s(le|LE)\s/},
    {token: "special", regex: /\s(ge|GE)\s/},
    {token: "special", regex: /\s(gt|GT)\s/},
    {token: "special", regex: /\s(or|OR)\s/},
    {token: "special", regex: /\s(and|AND)\s/},

    {token: "atom", regex: /"/},

    {token: "number", regex: /[-+]?(?:\.\d+|\d+\.?\d*)(?:e[-+]?\d+)?/i},

    {token: "operator", regex: /[-+\/*,\(\)]+/},

    {token: "variable", regex: /[a-zA-Z@:][\w@:$]*/},

  ]
});


//
// definition of the ".py" mode to be used by CodeMirror
//
CodeMirror.defineMode("python", function(conf, parserConf) {
    var ERRORCLASS = "error";

    function wordRegexp(words) {
        return new RegExp("^((" + words.join(")|(") + "))\\b");
    }

    var wordOperators = wordRegexp(["and", "or", "not", "is"]);
    var commonKeywords = ["as", "assert", "break", "class", "continue",
                          "def", "del", "elif", "else", "except", "finally",
                          "for", "from", "global", "if", "import",
                          "lambda", "pass", "raise", "return",
                          "try", "while", "with", "yield", "in"];
    var commonBuiltins = ["abs", "all", "any", "bin", "bool", "bytearray", "callable", "chr",
                          "classmethod", "compile", "complex", "delattr", "dict", "dir", "divmod",
                          "enumerate", "eval", "filter", "float", "format", "frozenset",
                          "getattr", "globals", "hasattr", "hash", "help", "hex", "id",
                          "input", "int", "isinstance", "issubclass", "iter", "len",
                          "list", "locals", "map", "max", "memoryview", "min", "next",
                          "object", "oct", "open", "ord", "pow", "property", "range",
                          "repr", "reversed", "round", "set", "setattr", "slice",
                          "sorted", "staticmethod", "str", "sum", "super", "tuple",
                          "type", "vars", "zip", "__import__", "NotImplemented",
                          "Ellipsis", "__debug__"];
    CodeMirror.registerHelper("hintWords", "python", commonKeywords.concat(commonBuiltins));

    function top(state) {
        return state.scopes[state.scopes.length - 1];
    }

    var delimiters = parserConf.delimiters || parserConf.singleDelimiters || /^[\(\)\[\]\{\}@,:`=;\.\\]/;
    //               (Backwards-compatiblity with old, cumbersome config system)
    var operators = [parserConf.singleOperators, parserConf.doubleOperators, parserConf.doubleDelimiters, parserConf.tripleDelimiters,
                     parserConf.operators || /^([-+*/%\/&|^]=?|[<>=]+|\/\/=?|\*\*=?|!=|[~!@])/]
    for (var i = 0; i < operators.length; i++) if (!operators[i]) operators.splice(i--, 1)

    var hangingIndent = parserConf.hangingIndent || conf.indentUnit;

    var myKeywords = commonKeywords, myBuiltins = commonBuiltins;
    if (parserConf.extra_keywords != undefined)
        myKeywords = myKeywords.concat(parserConf.extra_keywords);

    if (parserConf.extra_builtins != undefined)
        myBuiltins = myBuiltins.concat(parserConf.extra_builtins);

    var py3 = !(parserConf.version && Number(parserConf.version) < 3)
    if (py3) {
        // since http://legacy.python.org/dev/peps/pep-0465/ @ is also an operator
        var identifiers = parserConf.identifiers|| /^[_A-Za-z\u00A1-\uFFFF][_A-Za-z0-9\u00A1-\uFFFF]*/;
        myKeywords = myKeywords.concat(["nonlocal", "False", "True", "None", "async", "await"]);
        myBuiltins = myBuiltins.concat(["ascii", "bytes", "exec", "print"]);
        var stringPrefixes = new RegExp("^(([rbuf]|(br)|(fr))?('{3}|\"{3}|['\"]))", "i");
    } else {
        var identifiers = parserConf.identifiers|| /^[_A-Za-z][_A-Za-z0-9]*/;
        myKeywords = myKeywords.concat(["exec", "print"]);
        myBuiltins = myBuiltins.concat(["apply", "basestring", "buffer", "cmp", "coerce", "execfile",
                                        "file", "intern", "long", "raw_input", "reduce", "reload",
                                        "unichr", "unicode", "xrange", "False", "True", "None"]);
        var stringPrefixes = new RegExp("^(([rubf]|(ur)|(br))?('{3}|\"{3}|['\"]))", "i");
    }
    var keywords = wordRegexp(myKeywords);
    var builtins = wordRegexp(myBuiltins);

    // tokenizers
    function tokenBase(stream, state) {
        var sol = stream.sol() && state.lastToken != "\\"
        if (sol) state.indent = stream.indentation()
        // Handle scope changes
        if (sol && top(state).type == "py") {
            var scopeOffset = top(state).offset;
            if (stream.eatSpace()) {
                var lineOffset = stream.indentation();
                if (lineOffset > scopeOffset)
                    pushPyScope(state);
                else if (lineOffset < scopeOffset && dedent(stream, state) && stream.peek() != "#")
                    state.errorToken = true;
                return null;
            } else {
                var style = tokenBaseInner(stream, state);
                if (scopeOffset > 0 && dedent(stream, state))
                    style += " " + ERRORCLASS;
                return style;
            }
        }
        return tokenBaseInner(stream, state);
    }

    function tokenBaseInner(stream, state) {
        if (stream.eatSpace()) return null;

        // Handle Comments
        if (stream.match(/^#.*/)) return "comment";

        // Handle Number Literals
        if (stream.match(/^[0-9\.]/, false)) {
            var floatLiteral = false;
            // Floats
            if (stream.match(/^[\d_]*\.\d+(e[\+\-]?\d+)?/i)) { floatLiteral = true; }
            if (stream.match(/^[\d_]+\.\d*/)) { floatLiteral = true; }
            if (stream.match(/^\.\d+/)) { floatLiteral = true; }
            if (floatLiteral) {
                // Float literals may be "imaginary"
                stream.eat(/J/i);
                return "number";
            }
            // Integers
            var intLiteral = false;
            // Hex
            if (stream.match(/^0x[0-9a-f_]+/i)) intLiteral = true;
            // Binary
            if (stream.match(/^0b[01_]+/i)) intLiteral = true;
            // Octal
            if (stream.match(/^0o[0-7_]+/i)) intLiteral = true;
            // Decimal
            if (stream.match(/^[1-9][\d_]*(e[\+\-]?[\d_]+)?/)) {
                // Decimal literals may be "imaginary"
                stream.eat(/J/i);
                // TODO - Can you have imaginary longs?
                intLiteral = true;
            }
            // Zero by itself with no other piece of number.
            if (stream.match(/^0(?![\dx])/i)) intLiteral = true;
            if (intLiteral) {
                // Integer literals may be "long"
                stream.eat(/L/i);
                return "number";
            }
        }

        // Handle Strings
        if (stream.match(stringPrefixes)) {
            var isFmtString = stream.current().toLowerCase().indexOf('f') !== -1;
            if (!isFmtString) {
                // if this is a nested format string (e.g. f' {   f"{10*10}" + "a" }' )
                // we do not format the nested expression and treat the nested format
                // string as regular string
                state.tokenize = tokenStringFactory(stream.current());
                return state.tokenize(stream, state);
            } else {
                // need to do something more sophisticated
                state.tokenize = formatStringFactory(stream.current(), state.tokenize);
                return state.tokenize(stream, state);
            }
        }

        for (var i = 0; i < operators.length; i++)
            if (stream.match(operators[i])) return "operator"

        if (stream.match(delimiters)) return "punctuation";

        if (state.lastToken == "." && stream.match(identifiers))
            return "property";

        if (stream.match(keywords) || stream.match(wordOperators))
            return "keyword";

        if (stream.match(builtins))
            return "builtin";

        if (stream.match(/^(self|cls)\b/))
            return "variable-2";

        if (stream.match(identifiers)) {
            if (state.lastToken == "def" || state.lastToken == "class")
                return "def";
            return "variable";
        }

        // Handle non-detected items
        stream.next();
        return ERRORCLASS;
    }

    function formatStringFactory(delimiter, tokenOuter) {
        while ("rubf".indexOf(delimiter.charAt(0).toLowerCase()) >= 0)
            delimiter = delimiter.substr(1);

        var singleline = delimiter.length == 1;
        var OUTCLASS = "string";

        function tokenFString(stream, state) {
            // inside f-str Expression
            if (stream.match(delimiter)) {
                // expression ends pre-maturally, but very common in editing
                // Could show error to remind users to close brace here
                state.tokenize = tokenString
                return OUTCLASS;
            } else if (stream.match('{')) {
                // starting brace, if not eaten below
                return "punctuation";
            } else if (stream.match('}')) {
                // return to regular inside string state
                state.tokenize = tokenString
                return "punctuation";
            } else {
                // use tokenBaseInner to parse the expression
                return tokenBaseInner(stream, state);
            }
        }

        function tokenString(stream, state) {
            while (!stream.eol()) {
                stream.eatWhile(/[^'"\{\}\\]/);
                if (stream.eat("\\")) {
                    stream.next();
                    if (singleline && stream.eol())
                        return OUTCLASS;
                } else if (stream.match(delimiter)) {
                    state.tokenize = tokenOuter;
                    return OUTCLASS;
                } else if (stream.match('{{')) {
                    // ignore {{ in f-str
                    return OUTCLASS;
                } else if (stream.match('{', false)) {
                    // switch to nested mode
                    state.tokenize = tokenFString
                    if (stream.current()) {
                        return OUTCLASS;
                    } else {
                        // need to return something, so eat the starting {
                        stream.next();
                        return "punctuation";
                    }
                } else if (stream.match('}}')) {
                    return OUTCLASS;
                } else if (stream.match('}')) {
                    // single } in f-string is an error
                    return ERRORCLASS;
                } else {
                    stream.eat(/['"]/);
                }
            }
            if (singleline) {
                if (parserConf.singleLineStringErrors)
                    return ERRORCLASS;
                else
                    state.tokenize = tokenOuter;
            }
            return OUTCLASS;
        }
        tokenString.isString = true;
        return tokenString;
    }

    function tokenStringFactory(delimiter) {
        while ("rubf".indexOf(delimiter.charAt(0).toLowerCase()) >= 0)
            delimiter = delimiter.substr(1);

        var singleline = delimiter.length == 1;
        var OUTCLASS = "string";

        function tokenString(stream, state) {
            while (!stream.eol()) {
                stream.eatWhile(/[^'"\\]/);
                if (stream.eat("\\")) {
                    stream.next();
                    if (singleline && stream.eol())
                        return OUTCLASS;
                } else if (stream.match(delimiter)) {
                    state.tokenize = tokenBase;
                    return OUTCLASS;
                } else {
                    stream.eat(/['"]/);
                }
            }
            if (singleline) {
                if (parserConf.singleLineStringErrors)
                    return ERRORCLASS;
                else
                    state.tokenize = tokenBase;
            }
            return OUTCLASS;
        }
        tokenString.isString = true;
        return tokenString;
    }

    function pushPyScope(state) {
        while (top(state).type != "py") state.scopes.pop()
        state.scopes.push({offset: top(state).offset + conf.indentUnit,
                           type: "py",
                           align: null})
    }

    function pushBracketScope(stream, state, type) {
        var align = stream.match(/^([\s\[\{\(]|#.*)*$/, false) ? null : stream.column() + 1
        state.scopes.push({offset: state.indent + hangingIndent,
                           type: type,
                           align: align})
    }

    function dedent(stream, state) {
        var indented = stream.indentation();
        while (state.scopes.length > 1 && top(state).offset > indented) {
            if (top(state).type != "py") return true;
            state.scopes.pop();
        }
        return top(state).offset != indented;
    }

    function tokenLexer(stream, state) {
        if (stream.sol()) state.beginningOfLine = true;

        var style = state.tokenize(stream, state);
        var current = stream.current();

        // Handle decorators
        if (state.beginningOfLine && current == "@")
            return stream.match(identifiers, false) ? "meta" : py3 ? "operator" : ERRORCLASS;

        if (/\S/.test(current)) state.beginningOfLine = false;

        if ((style == "variable" || style == "builtin")
            && state.lastToken == "meta")
            style = "meta";

        // Handle scope changes.
        if (current == "pass" || current == "return")
            state.dedent += 1;

        if (current == "lambda") state.lambda = true;
        if (current == ":" && !state.lambda && top(state).type == "py")
            pushPyScope(state);

        if (current.length == 1 && !/string|comment/.test(style)) {
            var delimiter_index = "[({".indexOf(current);
            if (delimiter_index != -1)
                pushBracketScope(stream, state, "])}".slice(delimiter_index, delimiter_index+1));

            delimiter_index = "])}".indexOf(current);
            if (delimiter_index != -1) {
                if (top(state).type == current) state.indent = state.scopes.pop().offset - hangingIndent
                else return ERRORCLASS;
            }
        }
        if (state.dedent > 0 && stream.eol() && top(state).type == "py") {
            if (state.scopes.length > 1) state.scopes.pop();
            state.dedent -= 1;
        }

        return style;
    }

    var external = {
        startState: function(basecolumn) {
            return {
                tokenize: tokenBase,
                scopes: [{offset: basecolumn || 0, type: "py", align: null}],
                indent: basecolumn || 0,
                lastToken: null,
                lambda: false,
                dedent: 0
            };
        },

        token: function(stream, state) {
            var addErr = state.errorToken;
            if (addErr) state.errorToken = false;
            var style = tokenLexer(stream, state);

            if (style && style != "comment")
                state.lastToken = (style == "keyword" || style == "punctuation") ? stream.current() : style;
            if (style == "punctuation") style = null;

            if (stream.eol() && state.lambda)
                state.lambda = false;
            return addErr ? style + " " + ERRORCLASS : style;
        },

        indent: function(state, textAfter) {
            if (state.tokenize != tokenBase)
                return state.tokenize.isString ? CodeMirror.Pass : 0;

            var scope = top(state), closing = scope.type == textAfter.charAt(0)
            if (scope.align != null)
                return scope.align - (closing ? 1 : 0)
            else
                return scope.offset - (closing ? hangingIndent : 0)
        },

        electricInput: /^\s*[\}\]\)]$/,
        closeBrackets: {triples: "'\""},
        lineComment: "#",
        fold: "indent"
    };
    return external;
});
