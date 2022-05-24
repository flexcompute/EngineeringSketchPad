// ESP-gloves.js implements gloves functions for the Engineering Sketch Pad (ESP)
// written by Paul Mokotoff

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

// interface functions that ESP-gloves.js provides
//    gloves.launch()
//    gloves.cmdUndo()
//    gloves.cmdSolve()
//    gloves.cmdSave()
//    gloves.cmdQuit()
//
//    gloves.cmdHome()                      not provided
//    gloves.cmdLeft()                      not provided
//    gloves.cmdRite()                      not provided
//    gloves.cmdBotm()                      not provided
//    gloves.cmdTop()                       not provided
//    gloves.cmdIn()                        not provided
//    gloves.cmdOut()                       not provided
//
//    gloves.mouseDown(e)
//    gloves.mouseMove(e)
//    gloves.mouseUp(e)
//    gloves.mouseWheel(e)                  not provided
//    gloves.mouseLeftCanvas(e)
//
//    gloves.keyPress(e)
//    gloves.keyDown(e)                     not provided
//    gloves.keyUp(e)                       not provided
//    gloves.keyPressPart1(myKeyPress)
//    gloves.keyPressPart2(picking, gprim)
//    gloves.updateKeyWindow()
//
//    gloves.timLoadCB(text)
//    gloves.timSaveCB(text)
//    gloves.timQuitCB(text)
//    gloves.timMesgCB(text)
//
// additional functions:
//    gloves.addComponent(compType)
//    gloves.postMenu()
//    gloves.changeMode(newMode)
//    gloves.computeDist(t, m, xm, ym)
//    gloves.toggleCGs()
//

"use strict";

/* component categories for gloves.addDropdown */
gloves.PAYLOAD = 200;
gloves.AERO    = 300;
gloves.STRUCT  = 400;
gloves.CONTROL = 500;
gloves.PROP    = 600;
gloves.GEAR    = 700;
gloves.CUSTOM  = 800;

/* analysis categories for addAnalysis */
gloves.CG      = 1100;
gloves.AVL     = 1200;

/* curvature for BLENDing */
gloves.C2 = 4;
gloves.C1 = 2;
gloves.C0 = 1;


//
// callback when gloves is launched
//
gloves.launch = function () {
    // alert("in gloves.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // if there is no .csm file in session, ask user for a name
    var myFilename = null;
    if (wv.filenames == "|") {

        // wait until user gives a valid response
        while (myFilename === null) {
            myFilename = prompt("Enter the name of the .csm file that you want gloves to create");
        }

        // make sure it ends in ".csm"
        if (myFilename.length < 4 || myFilename.substring(myFilename.length-3) != ".csm") {
            myFilename += ".csm";
        }
    } else {
        var filenames = wv.filenames.split('|');
        myFilename = filenames[1];
    }

    // load the tim
    browserToServer("timLoad|gloves|"+myFilename+"|");

    // change done button legend and build the new menu
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Gloves";
    button.style.backgroundColor = "#3FFF3F";

    menu = document.getElementById("myDoneMenu");
    while (menu.firstChild != null) {
        menu.removeChild(menu.firstChild);
    }

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Toggle CGs";
    button.value   = "Toggle CGs";
    button.onclick = gloves.toggleCGs;
    menu.appendChild(button);

    menu.appendChild(document.createElement("hr"));

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Save and exit";
    button.value   = "Save";
    button.onclick = gloves.cmdSave;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Quit and exit";
    button.value   = "Quit";
    button.onclick = gloves.cmdQuit;
    menu.appendChild(button);

    // change solve button legend
    button = document.getElementById("solveButton");
    button["title"] = "Add Component to Configuration";
    button["innerHTML"] = "Add Component";
    button.style.backgroundColor = "#FFFF3F";

    // update the mode
    changeMode(9);
    gloves.changeMode(0);

    /* objects for storing design parameters and sensitivities */
    gloves.comp = [];

    /* remember the text canvas */
    gloves.cxt = document.getElementById("glovesText").getContext("2d");

    if (!gloves.cxt) {
        console.log("ERROR - no 2D canvas context.");
        return;
    }

    /* no components or pips selected yet */
    gloves.icomp = -1;
    gloves.ipip  = -1;
    gloves.iopt  = -1;

    gloves.moreUpd = 0;

    /* assume CG pips won't be drawn */
    gloves.drawCg = 0;

};


//
// callback when "gloves->Undo" is pressed (called by ESP.html)
//
gloves.cmdUndo = function () {
    alert("gloves.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
gloves.cmdSolve = function () {
    // alert("in gloves.cmdSolve()");

    if (gloves.mode == -1) {
        alert("ERROR - cannot add component during update.");
        return 1;

    }

    // determine which component to add
    var ans = prompt("What type of component to you want to add?\n"+
                     "Payload system:\n"+
                     "  11   Box\n"+
                     "  12   Cylinder\n"+
                     "  13   Sphere\n"+
                     "Aerodynamic system:\n"+
                     "  21   Fuselage\n"+
                     "  22   Wing\n"+
                     "  23   Vtail")
    if (ans === null) {
        return;
    } else if (ans == "11") {
        gloves.addComponent("glovesBox");
    } else if (ans == "12") {
        gloves.addComponent("glovesCyl");
    } else if (ans == "13") {
        gloves.addComponent("glovesSphr");
    } else if (ans == "21") {
        gloves.addComponent("glovesFuse");
    } else if (ans == "22") {
        gloves.addComponent("glovesWing");
    } else if (ans == "23") {
        gloves.addComponent("glovesVtail");
    } else {
        alert("Unknown component type");
    }
};


//
// callback when "gloves->Save" is pressed (called by ESP.html)
//
gloves.cmdSave = function () {
    // alert("in gloves.cmdSave()");

    if (gloves.mode == -1) {
        alert("ERROR - cannot save/exit GLOVES during update.");
        return 1;

    }

    /* clear the text canvas */
    var cxt = gloves.cxt;
    cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

    // return the done menu to its original state
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Done";
    button.style.backgroundColor = null;

    var menu = document.getElementById("myDoneMenu");
    while (menu.firstChild != null) {
        menu.removeChild(menu.firstChild);
    }

    // return the solve button to its original state
    button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    /* change solve button legend */
    var button = document.getElementById("solveButton");
    button["title"] = "Build or Solve";
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = "#FFFFFF";

    // execute save in the tim
    browserToServer("timSave|gloves|");

};


//
// callback when "gloves->Quit" is pressed (called by ESP.html)
//
gloves.cmdQuit = function () {
    // alert("in gloves.cmdQuit()");

    if (gloves.mode == -1) {
        alert("ERROR - cannot save/exit GLOVES during update.");
        return 1;

    }

    /* clear the text canvas */
    var cxt = gloves.cxt;
    cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

    // return the done menu to its original state
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Done";
    button.style.backgroundColor = null;

    var menu = document.getElementById("myDoneMenu");
    while (menu.firstChild != null) {
        menu.removeChild(menu.firstChild);
    }

    // return the solve button to its original state
    button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    /* change solve button legend */
    var button = document.getElementById("solveButton");
    button["title"] = "Build or Solve";
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = "#FFFFFF";

    // execute quit in the tim
    browserToServer("timQuit|gloves|");

};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//gloves.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//gloves.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//gloves.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//gloves.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//gloves.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//gloves.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//gloves.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
gloves.mouseDown = function(e) {

    if ((gloves.mode <= 0) ||
        (gloves.mode == 3) ) {
        main.mouseDown(e);

    } else {

        /* record mouse position to avoid accident translation */
        wv.startX =  e.clientX - wv.offLeft             - 1;
        wv.startY = -e.clientY + wv.offTop  + wv.height + 1;

    }
};


//
// callback when any mouse moves in canvas
//
gloves.mouseMove = function(e) {

    /* get the mouse position */
    gloves.mouseX = e.clientX - wv.offLeft;
    gloves.mouseY = e.clientY - wv.offTop ;

    if       ((gloves.mode <= 0) ||
              (gloves.mode == 3) ) {
        main.mouseMove(e);

    } else if (gloves.mode == 1) {

        /* get the menu */
        var menu = gloves.menu;

        /* assume mouse isn't over the menu */
        gloves.iopt = -1;

        /* find selected menu item */
        for (var item = 0; item < menu.length; item++) {
            if (gloves.mouseX >= menu[item].xmin-5 &&
                gloves.mouseX <= menu[item].xmax+5 &&
                gloves.mouseY >= menu[item].ymin-2 &&
                gloves.mouseY <= menu[item].ymax+2 ) {

                gloves.iopt = item;
                gloves.postMenu();
                break;

            }
        }

        if (gloves.iopt < 0) {

            /* clear the menu */
            var cxt = gloves.cxt;
            cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

            /* change the mode */
            gloves.changeMode(0);

            /* forget selected component */
            gloves.icomp = -1;
            gloves.ipip  = -1;

            /* forget menu and location */
            gloves.menu  = [];
            gloves.menuX = -1;
            gloves.menuY = -1;

            /* re-order the canvases */
            document.getElementById("WebViewer" ).style.zIndex = "1";
            document.getElementById("glovesText").style.zIndex = "0";

        }

    } else if (gloves.mode == 2) {

        /* get the cursor position */
        wv.cursorX =  e.clientX - wv.offLeft             - 1;
        wv.cursorY = -e.clientY + wv.offTop  + wv.height + 1;

        /* get the mouse screen coordinates */
        var xm = 2.0 * wv.cursorX / (wv.width  - 1.0) - 1.0;
        var ym = 2.0 * wv.cursorY / (wv.height - 1.0) - 1.0;

        /* get the transformation matrix */
        var m = wv.wireMatrix.getAsArray();

        /* setup the golden section search */
        var iter  = 0;

        var t0 = gloves.tt;
        var f0 = gloves.computeDist(t0, m, xm, ym);

        var tl = gloves.tt - gloves.delta;
        var tu = gloves.tt + gloves.delta;

        /* golden section search, part 1 */
        while (iter < 50) {
//          console.log("Part 1, iter = " + iter);
            /* evaluate candidate bounds */
            var fl = gloves.computeDist(tl, m, xm, ym);
            var fu = gloves.computeDist(tu, m, xm, ym);

            /* check if minimum was bounded */
            if ((fl - f0 > 1.0e-6) && (fu - f0 > 1.0e-6)) {
                break;

            } else {

                /* update bounds by a golden section */
                gloves.delta *= 1.618;

                tl = gloves.tt - gloves.delta;
                tu = gloves.tt + gloves.delta;

                /* iterate */
                iter++;

            }
        }

        console.log("Part 1 iterations = " + iter);

        /* compute the interval of uncertainty */
        var uncert = tu - tl;

        iter = 0;

        /* golden section search, part 2 */
        while ((uncert > 1.0e-6) && (iter < 100)) {
//          console.log("Part 2, iter = " + iter);

            /* get points in the interval */
            var ta = tl + 0.382 * uncert;
            var tb = tu - 0.382 * uncert;

            /* evaluate the points */
            var fa = gloves.computeDist(ta, m, xm, ym);
            var fb = gloves.computeDist(tb, m, xm, ym);

            /* reduce the interval of uncertainty */
            if (fa - fb > 1.0e-6) {
                tl = ta;
                fl = fa;

            } else {
                tu = tb;
                fu = fb;

            }

            /* compute the interval of uncertainty */
            uncert = tu - tl;

            /* iterate */
            iter++;

        }

        console.log("Part 2 iterations = " + iter);

        /* check that solution converged */
        if (iter == 100) {
            /* keep same topt */
            var topt = gloves.tt;

        } else {
            /* return the optimum t */
            var topt = (tl + tu) / 2;

        }
        console.log("topt = " + topt);

        /* update the value */
        var tempVal = gloves.curVal + topt;

        /* check if within bounds */
        if (tempVal > gloves.lbnd && tempVal < gloves.ubnd) {
            gloves.tt = topt;

        } else if (tempVal <= gloves.lbnd) {
            gloves.tt = gloves.lbnd - gloves.curVal;
            tempVal = gloves.lbnd;

        } else if (tempVal >= gloves.ubnd) {
            gloves.tt = gloves.ubnd - gloves.curVal;
            tempVal = gloves.ubnd;

        } else {
            gloves.tt = 0;

        }

        /* update rendering and pip */
        browserToServer("timMesg|gloves|newT|" + gloves.tt + "|");

        /* print value to screen */
        var cxt = gloves.cxt;
        cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

        gloves.tempVal = tempVal;

        var myText = gloves.curVar + " = " + gloves.tempVal;
        cxt.fillStyle = "black";
        cxt.fillText(myText, 6, 12);

    }
};


//
// callback when the mouse is released in canvas
//
gloves.mouseUp = function(e) {

    /* get the mouse position */
    gloves.mouseX = e.clientX - wv.offLeft;
    gloves.mouseY = e.clientY - wv.offTop ;

    if       ((gloves.mode <= 0) ||
              (gloves.mode == 3) ) {
        main.mouseUp(e);

    } else if (gloves.mode == 1) {

        /* get the menu */
        var menu = gloves.menu;

        /* find selected item */
        for (var item = 0; item < menu.length; item++) {
            if (gloves.mouseX >= menu[item].xmin-5 &&
                gloves.mouseX <= menu[item].xmax+5 &&
                gloves.mouseY >= menu[item].ymin-2 &&
                gloves.mouseY <= menu[item].ymax+2 ) {

                gloves.iopt = item;
                break;

            }
        }

        /* get context */
        var cxt = gloves.cxt;

        /* clear screen */
        cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

        /* keep the selected pip */
        var icomp = gloves.icomp;
        var ipip  = gloves.ipip ;
        var myPip = gloves.comp[icomp].pips[ipip];

        /* remember the branch */
        gloves.curBrch = gloves.comp[icomp].ibrch;

        /* remember the variable name */
        gloves.curVar = menu[item].labl;

        /* find parenthesis */
        var ibeg = gloves.curVar.search("\\(");

        if (ibeg == -1) {
            /* only one element in the parameter */
            gloves.varName = gloves.curVar;
            gloves.varIndx = 1;

        } else {
            /* multi-element parameter */
            var iend = gloves.curVar.search("\\)");

            /* get the name of the variable and its index */
            gloves.varName = gloves.curVar.substring(0     , ibeg);
            gloves.varIndx = gloves.curVar.substring(ibeg+1, iend);

        }

        /* check for cross-section operation */
        if ((gloves.curVar[0] == "+") ||
            (gloves.curVar[0] == "-") ||
            (gloves.curVar[0] == "*") ) {

            if         (gloves.varName.includes("Here")) {
                browserToServer("timMesg|gloves|sect|" + gloves.curBrch + "|" + gloves.varIndx + "|");

                /* change the mode */
                gloves.changeMode(-1);

            } else if  (gloves.varName.includes("NodeBefore")) {
                browserToServer("timMesg|gloves|node|" + gloves.curBrch + "|0|");

                /* change the mode */
                gloves.changeMode(-1);

            } else if ((gloves.varName.includes("SectBefore")) || (gloves.varName.includes("Before"))) {
                gloves.varIndx--;

                browserToServer("timMesg|gloves|sect|" + gloves.curBrch + "|" + gloves.varIndx + "|");

                /* change the mode */
                gloves.changeMode(-1);

            } else if  (gloves.varName.includes("NodeAfter")) {
                browserToServer("timMesg|gloves|node|" + gloves.curBrch + "|" + gloves.varIndx + "|");

                /* change the mode */
                gloves.changeMode(-1);

            } else if ((gloves.varName.includes("SectAfter")) || (gloves.varName.includes("After"))) {
                browserToServer("timMesg|gloves|sect|" + gloves.curBrch + "|" + gloves.varIndx + "|");

                /* change the mode */
                gloves.changeMode(-1);

            } else if  (gloves.varName.includes("Delete")) {
                browserToServer("timMesg|gloves|sect|" + gloves.curBrch + "|-" + gloves.varIndx + "|"); // section to be removed is gloves.varIndx-1

                /* change the mode */
                gloves.changeMode(-1);

            } else if  (gloves.varName.includes("Curvature")) {
                browserToServer("timMesg|gloves|curv|" + gloves.curBrch + "|" + gloves.varIndx + "|");

                /* change the mode */
                gloves.changeMode(-1);

            } else if  (gloves.varName.includes("Airfoil")) {

                /* get desired airfoil */
                var text = "Input the 4-digit NACA airfoil code (mptt):";
                var temp = Math.floor(Number(prompt(text, "0012")));

                if ((temp == null) || (isNaN(temp))) {
                    postMessage("Invalid airfoil selection: " + temp);

                } else if (temp > 9999) {
                    postMessage("Invalid airfoil selection (" + temp + ") - must be of the form: mptt");

                } else if (temp <= 0) {
                    postMessage("Invalid airfoil selection: " + temp);

                } else {

                    /* get the camber */
                    var m = Math.floor(temp / 1000);

                    /* "subtract" the camber */
                    temp -= 1000 * m;

                    /* get the maxloc */
                    var p = Math.floor(temp / 100);

                    /* "subtract" the maxloc */
                    temp -= 100 * p;

                    /* get the thickness */
                    var tt = temp / 100;

                    /* update parameters */
                    browserToServer("setPmtr|" + gloves.comp[icomp].name + ":camber" + "|" + gloves.varIndx + "|1|" +  m/100 + "|");
                    browserToServer("setPmtr|" + gloves.comp[icomp].name + ":maxloc" + "|" + gloves.varIndx + "|1|" +  p/10  + "|");
                    browserToServer("setPmtr|" + gloves.comp[icomp].name + ":thick"  + "|" + gloves.varIndx + "|1|" + tt     + "|");

                    /* change the airfoil */
                    browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

                    /* change the mode */
                    gloves.changeMode(-1);

                }

            } else if  (gloves.varName.includes("ChngHalf")) {

                /* get desired cylinder half (or whole) */
                var text  = "Input the desired cylinder half:\n";
                text += "    0 - both halves\n";
                text += "    1 - lower (-z) half\n";
                text += "    2 - upper (+z) half\n";

                var temp = Math.floor(Number(prompt(text, "0")));

                if ((temp == null) || (isNaN(temp))) {
                    postMessage("Invalid cylinder half: " + temp);

                } else if (temp > 2) {
                    postMessage("Invalid cylinder half (" + temp + ") - must be 0, 1, or 2");

                } else if (temp < 0) {
                    postMessage("Invalid cylinder half (" + temp + ") - must be 0, 1, or 2");

                } else {

                    /* change the cylinder half */
                    browserToServer("setPmtr|" + gloves.comp[icomp].name + ":half|1|1|" + temp + "|");
                    browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

                    /* change the mode */
                    gloves.changeMode(-1);

                }
            }

        } else {

            /* remember the body name */
            var compName = menu[item].labl.split(":")[0];
            gloves.compName = compName;

            /* manipulate the component */
            browserToServer("timMesg|gloves|chngComp|" + gloves.comp[icomp].ibrch + "|" + gloves.varName + "|" + gloves.varIndx + "|" + myPip.entIndx + "|" + myPip.entTess + "|");

            var myText = gloves.curVar + " = ";
            cxt.fillStyle = "black";
            cxt.fillText(myText, 6, 12);

            /* re-order the canvases */
            var webCanvas = document.getElementById("WebViewer");
            var txtCanvas = document.getElementById("glovesText");

            webCanvas.style.zIndex = "1";
            txtCanvas.style.zIndex = "0";

            /* change the mode */
            gloves.changeMode(2);

            /* turn off all pips except selected one */
            for (var jpip = 1; jpip <= gloves.comp[icomp].pips.length; jpip++) {

                if (jpip == ipip+1) continue;

                /* get the pip name */
                var pipName = compName + " Pip " + jpip;

                /* turn off the pip */
                wv.sceneGraph[pipName].attrs &= ~wv.plotAttrs.ON;

            }

            /* change the solve button */
            var button = document.getElementById("solveButton");
            button["innerHTML"] = "Updating ...";

            /* initial "t" is 0 (current location) */
            gloves.tt = 0;

            /* initial "delta" for golden section search */
            gloves.delta = 0.1;

        }

        /* forget the menu and its location */
        gloves.menu  = [];
        gloves.menuX = -1;
        gloves.menuY = -1;

    } else if (gloves.mode == 2) {

        /* clear the text */
        var cxt = gloves.cxt;

        cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

        /* turn on all pips */
        for (var jpip = 1; jpip <= gloves.comp[gloves.icomp].pips.length; jpip++) {

            if (jpip == ipip+1) continue;

            /* get the pip name */
            var pipName = gloves.compName + " Pip " + jpip;

            /* turn on the pip */
            wv.sceneGraph[pipName].attrs |= wv.plotAttrs.ON;

        }

        /* find the design parameter index */
        var myKey = gloves.curVar;

        /* update the component */
        browserToServer("setPmtr|" + gloves.varName + "|" + gloves.varIndx + "|1|" + gloves.tempVal + "|");
        browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

        /* forget pip location and sensitivity */
        gloves.xx = gloves.yy = gloves.zz = 0;
        gloves.dx = gloves.dy = gloves.dz = 0;

        /* forget bounds */
        gloves.lbnd = gloves.ubnd = 0;

        /* forget search parameters */
        gloves.tt    = 0  ;
        gloves.delta = 0.1;

        /* forget parameter value */
        gloves.curVal = gloves.tempVal = 0;
        gloves.varIndex = 0;

        /* forget variable name */
        gloves.varName = null;
        gloves.curVar  = null;

        /* change the mode */
        gloves.changeMode(-1);

    } else if (gloves.mode == 3) {

        console.log("in gloves.mode = 3");

    }
};


//
// callback when the mouse wheel is rolled in canvas
//
//gloves.mouseWheel = function(e) {
//    main.mouseWheel(e);
//};


//
// callback when the mouse leaves the canvas
//
gloves.mouseLeftCanvas = function(e) {

    /* post a message */
    if (gloves.mode == 2) {
        postMessage("Mouse left the canvas --> stopped updating component at" + gloves.curVar + " = " + gloves.tempVal);
        gloves.mouseUp(e);
    }

};


//
// callback when a key is pressed
//
gloves.keyPress = function (e) {
    // alert("in gloves.keyPress");

    /* quit updating the component */
    if        (e.key == "q" && gloves.mode == 2) {

        /* return BODY to original configuration */
        browserToServer("timMesg|gloves|newT|0|");

        /* change the mode */
        gloves.changeMode(0);

        /* clear the text */
        var cxt = gloves.cxt;

        cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

        /* turn on all pips */
        var ipip = gloves.ipip;

        for (var jpip = 1; jpip <= gloves.comp[gloves.icomp].pips.length; jpip++) {

            if (jpip == ipip+1) continue;

            /* get the pip name */
            var pipName = gloves.compName + " Pip " + jpip;

            /* turn on the pip */
            wv.sceneGraph[pipName].attrs |= wv.plotAttrs.ON;

        }

        /* forget pip location and sensitivity */
        gloves.xx = gloves.yy = gloves.zz = 0;
        gloves.dx = gloves.dy = gloves.dz = 0;

        /* forget bounds */
        gloves.lbnd = gloves.ubnd = 0;

        /* forget search parameters */
        gloves.tt    = 0  ;
        gloves.delta = 0.1;

        /* forget selected component */
        gloves.icomp   = gloves.ipip    = gloves.iopt = -1;
        gloves.curBrch =                                    -1;
        gloves.curVal  = gloves.tempVal =                  0;

        /* forget variable name */
        gloves.curVar = null;

        /* change the mode */
        gloves.changeMode(0);

        return 1;

    }

    /* update the value of the parameter */
    if        (e.key == "v" && gloves.mode == 2) {

        /* prompt user for a value */
        var text = "Input value for " + gloves.curVar + ":";
        var temp = prompt(text, gloves.curVal);

        if (temp != null && temp != NaN) {

            /* convert input to a number */
            temp = Number(temp);

            if (temp > gloves.lbnd && temp < gloves.ubnd) {

                /* find the design parameter index */
                var myKey = gloves.curVar;

                /* update the component */
                browserToServer("setPmtr|" + gloves.varName + "|" + gloves.varIndx + "|1|" + temp + "|");
                browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

                /* change the mode */
                gloves.changeMode(-1);

            } else {

                postMessage("Returning to configuration before adjusting component.")

                /* return BODY to original configuration */
                browserToServer("timMesg|gloves|newT|0|");

                /* change the mode */
                gloves.changeMode(0);

            }

        } else {

            postMessage("Returning to configuration before adjusting component.")

            /* return BODY to original configuration */
            browserToServer("timMesg|gloves|newT|0|");

            /* change the mode */
            gloves.changeMode(0);

        }

        /* clear the text */
        var cxt = gloves.cxt;

        cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

        /* turn on all pips */
        var ipip = gloves.ipip;

        for (var jpip = 1; jpip <= gloves.comp[gloves.icomp].pips.length; jpip++) {

            if (jpip == ipip+1) continue;

            /* get the pip name */
            var pipName = gloves.compName + " Pip " + jpip;

            /* turn on the pip */
            wv.sceneGraph[pipName].attrs |= wv.plotAttrs.ON;

        }

        /* forget pip location and sensitivity */
        gloves.xx = gloves.yy = gloves.zz = 0;
        gloves.dx = gloves.dy = gloves.dz = 0;

        /* forget the bounds */
        gloves.lbnd = gloves.ubnd = 0;

        /* forget parameter value */
        gloves.curVal = gloves.tempVal = 0;

        /* forget search parameters */
        gloves.tt    = 0  ;
        gloves.delta = 0.1;

        /* forget variable name */
        gloves.curVar = null;

        return 1;

    }

    /* complete parent-child connections */
    if ((e.key == "d") && (gloves.mode == 3)) {

        /* return to original mode */
        gloves.changeMode(0);

        if (gloves.ichild.length > 0) {

            /* create message */
            var mesg = "timMesg|gloves|linkComp|" + gloves.curBrch + "|";

            for (var irow = 0; irow < gloves.ichild.length; irow++) {
                mesg += gloves.ichild[irow] + "|";

            }

            /* send message to the browser */
            browserToServer(mesg);

        }

        /* free array */
        gloves.iparent =  0 ;
        gloves.ichild  = [ ];

        return 1;

    }

    /* if key press is not handled */
    return 0;

};


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//gloves.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//gloves.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
gloves.keyPressPart1 = function(myKeyPress) {

    // alert("in gloves.keyPressPart1(myKeyPress="+myKeyPress+")");

    /* pop-up a menu for the pip */
    if ((myKeyPress == "m") && (gloves.mode == 0)) {
        wv.pick    =  1;
        wv.picking = 77;

        return 1;
    }

    /* create parent-child connections */
    if ((myKeyPress == "c") && (gloves.mode == 0)) {
        wv.pick    =  1;
        wv.picking = 67;

        return 1;
    }

    /* add parent-child connection */
    if ((myKeyPress == "c") && (gloves.mode == 3)) {
        wv.pick    =  1;
        wv.picking = 67;

        return 1;
    }

    if ((myKeyPress == "r") && (gloves.mode == 3)) {
        wv.pick    =  1;
        wv.picking = 82;

        return 1;
    }

    /* key press wasn't handled */
    return 0;

};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
gloves.keyPressPart2 = function(picking, gprim) {
    // alert("in gloves.keyPressPart2(picking="+picking+"   gprim="+gprim+")");

    /* remember the mouse position */
    gloves.menuX = gloves.mouseX;
    gloves.menuY = gloves.mouseY;

    /* "c" - create parent-child connection */
    if ((picking == 67) && (gloves.mode == 3)) {

        /* get the component name */
        var compName = wv.picked.gprim.split(" ")[0];

        /* find the component index */
        for (var icomp = 0; icomp < gloves.comp.length; icomp++) {
            if (compName == gloves.comp[icomp].name) break;
        }

        if (icomp == gloves.comp.length) {
            postMessage("ERROR - component not found ('c')");
            return;
        }

        /* get length of children array */
        var nrow = gloves.ichild.length;

        /* check if component is already marked */
        for (var irow = 0; irow < nrow; irow++) {
            if (gloves.comp[icomp].ibrch == gloves.ichild[irow]) break;
        }

        if (gloves.comp[icomp].ibrch == gloves.curBrch) {
            postMessage(gloves.comp[icomp].name + " can't be selected --- it is already the parent.");

            return;
        }

        /* if loop reaches nrow, add the child */
        if (irow == nrow) {
            gloves.ichild[nrow] = gloves.comp[icomp].name;

            /* write message */
            postMessage("Adding " + gloves.comp[icomp].name + " as a child");

        } else {
            postMessage(gloves.comp[icomp].name + " is already a child");
            return;

        }
    }

    /* "c" - initiate parent-child connections */
    if ((picking == 67) && (gloves.mode == 0)) {

        /* get the component name (separated by a space) */
        var compName = wv.picked.gprim.split(" ")[0];

        /* find the component index */
        for (var icomp = 0; icomp < gloves.comp.length; icomp++) {
            if (compName == gloves.comp[icomp].name) break;
        }

        if (icomp == gloves.comp.length) {
            postMessage("ERROR - component not found ('c')");
            return;
        }

        /* remember parent index */
        gloves.curBrch = gloves.comp[icomp].ibrch;

        /* create array for remembering branch indices */
        gloves.ichild = [];

        /* change the mode */
        gloves.changeMode(3);

    }

    /* "r" - remove parent-child connection */
    if (picking == 82) {

        /* get the component name */
        var compName = wv.picked.gprim.split(" ")[0];

        /* find the component index */
        for (var icomp = 0; icomp < gloves.comp.length; icomp++) {
            if (compName == gloves.comp[icomp].name) break;
        }

        if (icomp == gloves.comp.length) {
            postMessage("ERROR - component not found ('r')");
            return;
        }

        /* get length of children array */
        var nrow = gloves.ichild.length;

        /* check if component is already marked */
        for (var irow = 0; irow < nrow; irow++) {
            if (gloves.comp[icomp].ibrch == gloves.ichild[irow]) break;
        }

        /* if loop reaches nrow, the child wasn't connected */
        if (irow == nrow) {
            postMessage(gloves.comp[icomp].name + " wasn't a child");

        } else {
            gloves.ichild.splice(irow, 1);
            console.log(gloves.ichild);

            /* write a message */
            postMessage("Removing " + gloves.comp[icomp].name + " as a child");
        }
    }

    /* "m" - list the menu and change modes */
    if (picking == 77) {

        /* get the component name (separated by a space) */
        var compName = wv.picked.gprim.split(" ")[0];

        /* match the selected component */
        for (var icomp = 0; icomp < gloves.comp.length; icomp++) {

            /* get the first pip name */
            var pipName = gloves.comp[icomp].pips[0].pipName.split(" ")[0];

            if (pipName == compName) {
                break;
            }
        }

        /* return if pip wasn't selected */
        if (icomp == gloves.comp.length) return;

        /* remember the component index */
        gloves.icomp = icomp;

        /* find the pip on that component */
        for (var ipip = 0; ipip < gloves.comp[icomp].pips.length; ipip++) {
            if (wv.picked.gprim == gloves.comp[icomp].pips[ipip].pipName) {
                break;
            }
        }

        /* ignore remainder if pip wasn't selected */
        if (ipip == gloves.comp[icomp].pips.length) {
            gloves.icomp = -1;
            return;

        }

        /* remember the pip index */
        gloves.ipip = ipip;

        /* assume first option is selected */
        gloves.iopt = 0;

        /* post the menu */
        gloves.postMenu();

        /* change the mode */
        gloves.changeMode(1);

        /* place text canvas in-front of wv canvas */
        document.getElementById("WebViewer" ).style.zIndex = "0";
        document.getElementById("glovesText").style.zIndex = "1";

    }
};


//
// function to update the key window
//
gloves.updateKeyWindow = function () {

    var mesg = "gloves - Available Key Presses:\n";

    if        (gloves.mode == 0) {
        mesg += "  'c' --> Connect parent/child\n";
        mesg += "          components\n";
        mesg += "  'm' --> Open a pip menu\n";

    } else if (gloves.mode == 1) {
        mesg += "  None\n";

    } else if (gloves.mode == 2) {
        mesg += "  'q' --> Quit without changing\n";
        mesg += "  'v' --> Input value for parameter\n";

    } else if (gloves.mode == 3) {
        mesg += "  'c' --> Connect component to parent\n";
        mesg += "  'd' --> Done connecting\n";
        mesg += "  'r' --> Remove component as child\n";

    } else if (gloves.mode <  0) {
        mesg += "  Updating ... none available\n";

    } else {
        alert("ERROR - mode " + newMode + " is not supported.");
        return -1;

    }

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    var timStatus = document.getElementById("timStatus");

    timStatus.replaceChild(pre, timStatus.lastChild);

};


//
// callback when timLoad returns
//
gloves.timLoadCB = function (text) {

    /* if component(s) exist, fill structure */
    if (text.length > 15) {

        /* change the mode */
        gloves.changeMode(-1);

        /* parse the string */
        var nameBrch = text.substring(9, text.length-1).split("|");
        console.log(nameBrch);

        /* get the names and branches as comma-separated list */
        var name = nameBrch[0].substring(1, nameBrch[0].length-1);
        var brch = nameBrch[1].substring(1, nameBrch[1].length-1);

        /* split comma-separated list into individual entities */
        var compName = name.split(",");
        var compBrch = brch.split(",");

        var ncomp = compName.length;

        /* remember these entities */
        for (var icomp = 0; icomp < ncomp; icomp++) {

            /* add component to the structure */
            gloves.comp[icomp] = {};

            /* remember the branch index and name */
            gloves.comp[icomp].name  = compName[icomp];
            gloves.comp[icomp].ibrch = compBrch[icomp];

        }

        /* update first component (triggers others in updateEnd callback) */
        gloves.curBrch = gloves.comp[0].ibrch;

        if (ncomp == 1) {
            gloves.moreUpd = 0;

        } else {
            gloves.moreUpd = 1;

        }

        browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

    }
};


//
// callback when timSave returns
//
gloves.timSaveCB = function (text) {

    /* get the updated parameter/branch lists */
    browserToServer("getPmtrs|");
    browserToServer("getBrchs|");

    /* update the tree window */
    console.log("rebuildTree --> Gloves --> timSaveCB");
    rebuildTreeWindow();

    /* return to "normal" ESP mode */
    changeMode(0);

    /* free the components */
    gloves.comp = [];

    /* no gloves mode exists */
    gloves.mode = null;

};


//
// callback when timQuit returns
//
gloves.timQuitCB = function (text) {

    /* get the updated parameter/branch lists */
    browserToServer("getPmtrs|");
    browserToServer("getBrchs|");

    /* update the tree window */
    console.log("rebuildTree --> Gloves --> timQuitCB");
    rebuildTreeWindow();

    /* return to "normal" ESP mode */
    changeMode(0);

    /* free the components */
    gloves.comp = [];

    /* no gloves mode exists */
    gloves.mode = null;

};


//
// callback when timMesg returns
//
gloves.timMesgCB = function (text) {

    /* remove "gloves|" */
    text = text.substring(7);

    /* find the message */
    if        (text.substring(0, 8) == "addComp|") {

        /* add a component to the structure */
        var icomp = gloves.comp.length;
        gloves.comp[icomp]      = {};

        /* get arguments of string */
        var ibrch = Number(text.substring(8, text.length-1));

        /* remember its branch index */
        gloves.comp[icomp].ibrch = ibrch;
        gloves.curBrch           = ibrch;

        /* remember its name */
        gloves.comp[icomp].name  = gloves.compName;

        /* get the sensitivities */
        browserToServer("timMesg|gloves|updateComp|" + gloves.comp[icomp].ibrch + "|");

        /* change the mode */
        gloves.changeMode(-1);

    } else if (text.substring(0, 10) == "updateBeg|") {

        /* initialize the message */
        gloves.message = text.substring(10, text.length-1);

        /* request next part of message */
        browserToServer("timMesg|gloves|updateMid|");

    } else if (text.substring(0, 10) == "updateMid|") {

        /* add on to the message */
        gloves.message += text.substring(10, text.length-1);

        /* request next part of message */
        browserToServer("timMesg|gloves|updateMid|");

    } else if (text.substring(0, 10) == "updateEnd|") {

        /* find the appropriate branch */
        for (var icomp = 0; icomp < gloves.comp.length; icomp++) {
            if (gloves.comp[icomp].ibrch == gloves.curBrch) {
                break;
            }
        }

        /* parse the JSON string */
        gloves.comp[icomp].pips = JSON.parse(gloves.message);

        /* remove the message */
        gloves.message = null;

        /* components or pips no longer selected */
        gloves.icomp   = -1;
        gloves.ipip    = -1;
        gloves.iopt    = -1;
        gloves.curBrch = -1;

        if (gloves.moreUpd > 0) {

            /* get current branch */
            gloves.curBrch = gloves.comp[gloves.moreUpd].ibrch;
            console.log("Current Branch = " + gloves.curBrch);

            /* update next component */
            browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

            /* check for more components */
            if (gloves.moreUpd == gloves.comp.length-1) {
                gloves.moreUpd = 0;

            } else {
                gloves.moreUpd += 1;

            }

            console.log(gloves.moreUpd);

        } else {

            /* get the branches and parameters */
            browserToServer("getPmtrs|");
            browserToServer("getBrchs|");

            /* draw the components */
            browserToServer("timDraw|gloves|");
            browserToServer("timMesg|gloves|drawPips|" + gloves.drawCg + "|");

        }

    } else if (text.substring(0, 9) == "chngComp|") {

        /* intialize the message */
        gloves.message = text.substring(9, text.length-1);

        /* parse the message */
        var splitString = gloves.message.split("|");

        /* remember the value of the design parameter */
        gloves.curVal = Number(splitString[0]);

        /* remember the bounds */
        var sepVals = splitString[1].substring(1, splitString[1].length-1).split(",");

        gloves.lbnd = Number(sepVals[0]);
        gloves.ubnd = Number(sepVals[1]);

        /* split up the coordinates at the pip and remember them */
        var sepVals = splitString[2].substring(1, splitString[2].length-1).split(",");

        gloves.xx = Number(sepVals[0]);
        gloves.yy = Number(sepVals[1]);
        gloves.zz = Number(sepVals[2]);

        /* split up the sensitivity at the pip and remember them */
        var sepVals = splitString[3].substring(1, splitString[3].length-1).split(",");

        gloves.dx = Number(sepVals[0]);
        gloves.dy = Number(sepVals[1]);
        gloves.dz = Number(sepVals[2]);

    } else if (text.substring(0, 9) == "drawPips|") {

        /* change the mode */
        gloves.changeMode(0);

//      /* get the branches and parameters */
//      browserToServer("getPmtrs|");
//      browserToServer("getBrchs|");

        console.log(wv.sceneGraph);

        /* loop through all gprims to find FACEs */
        var nnode = myTree.name.length;

        for (var inode = 0; inode < nnode; inode++) {

            /* check the name */
            if ((myTree.gprim[inode] != "") && (myTree.gprim[inode] != null)) {
                if (myTree.gprim[inode].includes("Face")) {

                    console.log(myTree.gprim[inode]);

                    /* toggle the transparency on/off */
                    if (gloves.drawCg == 1) {
//                      wv.sceneGraph[myTree.gprim[inode]].attrs |=  wv.plotAttrs.TRANSPARENT;
                        changeProp(inode, 3, "on" );

                    } else {
//                      wv.sceneGraph[myTree.gprim[inode]].attrs &= ~wv.plotAttrs.TRANSPARENT;
                        changeProp(inode, 3, "off");

                    }
                }
            }
        }

        /* update the tree */
        myTree.update();

        console.log(wv.sceneGraph);

    } else if (text.substring(0, 5) == "newT|") {

        /* do nothing */

    } else if (text.substring(0, 5) == "sect|") {

        browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

        /* change the mode */
        gloves.changeMode(-1);

    } else if (text.substring(0, 5) == "node|") {

        browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

        /* change the mode */
        gloves.changeMode(-1);

    } else if (text.substring(0, 5) == "curv|") {

        /* read in the arguments */
        gloves.message = text.substring(5, text.length-1);
        console.log(gloves.message);

        /* parse the message */
        var splitString = gloves.message.split("|");
        console.log(splitString);

        /* remember the component name */
        gloves.varName = splitString[0];

        /* remember the component index */
        gloves.varIndx = splitString[1];

        /* remember the available options */
        gloves.curVar  = splitString[2];

        /* create message for prompt */
        var text = "Input integer from the available curvatures:\n";
        var errs = "ERROR: invalid curvature selected (must be ";

        /* check for valid curvatures */
        if (gloves.curVar & gloves.C2) {
            text += "    2 - C2 (curvature continuous)\n";
            errs += "2, ";

        }

        if (gloves.curVar & gloves.C1) {
            text += "    1 - C1 (slope continuous)\n";
            errs += "1, ";

        }

        if (gloves.curVar & gloves.C0) {
            text += "    0 - C0 (value continuous)\n";
            errs += "0";
        }

        errs += ")";

        var newCurv = parseInt(prompt(text, 2));
        console.log(newCurv);

        if        (isNaN(newCurv)) {
            postMessage("WARNING: new curvature was not selected.");
            return;

        } else if ((newCurv < 0) || (newCurv > 2)) {
            postMessage("ERROR: invalid curvature selected (must be between 0 and 2).");
            return;

        } else if (((gloves.curVar & gloves.C2) == 0) && (newCurv == 2)) {
            postMessage(errs);
            return;

        } else if (((gloves.curVar & gloves.C1) == 0) && (newCurv == 1)) {
            postMessage(errs);
            return;

        } else if (((gloves.curVar & gloves.C0) == 0) && (newCurv == 0)) {
            postMessage(errs);
            return;

        }

        /* update component */
        browserToServer("setPmtr|" + gloves.varName + "|" + gloves.varIndx + "|1|" + newCurv + "|");
        browserToServer("timMesg|gloves|updateComp|" + gloves.curBrch + "|");

        /* change the mode */
        gloves.changeMode(-1);

    } else if (text.substring(0, 9) == "linkComp|") {

        /* draw the components */
        browserToServer("timDraw|gloves|");
        browserToServer("timMesg|gloves|drawPips|" + gloves.drawCg + "|");

    } else {

        /* post to the message window */
        postMessage(text);

    }
};


//
// callback when user asks to add a component
//
gloves.addComponent = function(compType) {
    // alert("in gloves.addComponent(compType="+compType+")");

    /* request the name of the component */
    var compName = prompt("Input the name of the component: ", "tempComp" + gloves.comp.length);

    if (!compName) {
        return -1;
    }

    /* eliminate any excess spaces before/after name */
    compName = compName.replace(/ /g, "");

    /* request the component to be built */
    browserToServer("timMesg|gloves|addComp|" + compType + "|" + compName + "|");

    /* remember the name */
    gloves.compName = compName;

};


//
// callback to post menu when user presses "m" on a Pip
//
gloves.postMenu = function() {
    // alert("in gloves.postMenu()");

    /* get the menu */
    var menuList = gloves.comp[gloves.icomp].pips[gloves.ipip].pmtrs;

    /* clear the canvas */
    var cxt = gloves.cxt;
    cxt.clearRect(0, 0, cxt.canvas.width, cxt.canvas.height);

    /* place the menu */
    var xpos = gloves.menuX;
    var ypos = gloves.menuY;

    /* find equired menu width */
    var width = 0;

    for (var item = 0; item < menuList.length; item++) {
        var foo = cxt.measureText(menuList[item]).width;

        if (foo > width) {
            width = foo;
        }
    }

    /* menu height */
    var height = 16;

    /* make and post the menu */
    gloves.menu = [];

    for (item = 0; item < menuList.length; item++) {

        /* make menu item object */
        gloves.menu[item]      = {};
        gloves.menu[item].labl = menuList[item];
        gloves.menu[item].xmin = xpos;
        gloves.menu[item].xmax = xpos + width;
        gloves.menu[item].ymin = ypos + (item - 0.5) * height;
        gloves.menu[item].ymax = ypos + (item + 0.5) * height;

        /* post the menu */
        if (item == gloves.iopt) {

            cxt.fillStyle = "blue";
            cxt.fillRect(gloves.menu[item].xmin-2,
                         gloves.menu[item].ymin-1,
                         gloves.menu[item].xmax-gloves.menu[item].xmin+4,
                         gloves.menu[item].ymax-gloves.menu[item].ymin+2);

            cxt.fillStyle = "white";
            cxt.fillText(gloves.menu[item].labl,
                         gloves.menu[item].xmin,
                         gloves.menu[item].ymin+height/2);

        } else {

            cxt.fillStyle = "white";
            cxt.fillRect(gloves.menu[item].xmin-2,
                         gloves.menu[item].ymin-1,
                         gloves.menu[item].xmax-gloves.menu[item].xmin+4,
                         gloves.menu[item].ymax-gloves.menu[item].ymin+2);

            cxt.fillStyle = "black";
            cxt.fillText(gloves.menu[item].labl,
                         gloves.menu[item].xmin,
                         gloves.menu[item].ymin+height/2);

        }
    }

};


//
// function to change Mode within gloves
//
gloves.changeMode = function(newMode) {

    /* change the mode */
    gloves.mode = newMode;

    /* update the key window */
    gloves.updateKeyWindow();

    /* update the solve button */
    var button;

    if        (newMode == -1) {
        button = document.getElementById("solveButton");
        button["innerHTML"] = "Updating...";

    } else if (newMode ==  0) {
        button = document.getElementById("solveButton");
        button["innerHTML"] = "Add Component";

    }

    return 0;
};


//
// function to compute distance bwteeen points
//
gloves.computeDist = function(t, m, xm, ym) {

    /* get the pip coordinates */
    var x3d = gloves.xx + gloves.dx * t;
    var y3d = gloves.yy + gloves.dy * t;
    var z3d = gloves.zz + gloves.dz * t;

    /* scale by the focus */
    var xp = (x3d - wv.focus[0]) / wv.focus[3];
    var yp = (y3d - wv.focus[1]) / wv.focus[3];
    var zp = (z3d - wv.focus[2]) / wv.focus[3];

    /* convert to two-dimensional coordinates */
    var xtemp = m[0] * xp + m[4] * yp + m[ 8] * zp + m[12];
    var ytemp = m[1] * xp + m[5] * yp + m[ 9] * zp + m[13];
    var wtemp = m[3] * xp + m[7] * yp + m[11] * zp + m[15];

    /* convert to screen coordinates and compute the velocities */
    var x2d = xtemp / wtemp ;
    var y2d = ytemp / wtemp ;

    /* compute the square of the distance between the mouse and pip */
    var dx = x2d - xm;
    var dy = y2d - ym;

    var dist = dx * dx + dy * dy;

    return dist;

};


//
// callback to toggle the display of CGs
//
gloves.toggleCGs = function() {
    /* toggle CGs to be drawn */
    gloves.drawCg = 1 - gloves.drawCg;
    console.log("drawCg = " + gloves.drawCg);

    /* make call to server */
    browserToServer("timDraw|gloves|");
    browserToServer("timMesg|gloves|drawPips|" + gloves.drawCg + "|");

    /* change the mode */
    gloves.changeMode(-1);
};
