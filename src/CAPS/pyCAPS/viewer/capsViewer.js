// capsViewer.js designed to mimic ESP.js


// Post a message to the bottom frame 
function postMessage(mesg) {
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
}

function postParseError(text) {
    postMessage(" Server Message: Error parsing " + text);
}

// Send a message to the server
function browserToServer(text) {
    // alert("browserToServer(text="+text+")");

    if (wv.debugUI) {
        var date = new Date;
        console.log("("+date.toTimeString().substring(0,8)+") browser-->server: "+text.substring(0,40));
    }

    wv.socketUt.send(text);
}

//
// callback when "onresize" event occurs (called by capsViwer.html)
//
// resize the frames (with special handling to width of tlframe and height of brframe)
//
function resizeFrames() {
    // alert("resizeFrames()");

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
}


//---------------------------------------------------------------------------//
//
// Tree
//
//---------------------------------------------------------------------------//

function TreeAddNode(iparent, name, gprim, prop1, cbck1, prop2, cbck2, prop3, cbck3) {
    // alert("in TreeAddNode(" + iparent + "," + name + "," + gprim + "," + prop1 + "," + prop2 + "," + prop3 + ")");

    // validate the input
    if (iparent < 0 || iparent >= this.name.length) {
        alert("iparent=" + iparent + " is out of range");
        return;
    }

    // find the next node index
    var inode = this.name.length;

    // store this node's values
    this.name[  inode] = name;
    this.gprim[ inode] = gprim;
    this.parent[inode] = iparent;
    this.child[ inode] = -1;
    this.next[  inode] = -1;
    this.opened[inode] =  0;

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
    //    new node to the parent
    if (this.child[iparent] < 0) {
        this.child[iparent] = inode;

    // otherwise link this node to the last parent's child
    } else {
        var jnode = this.child[iparent];
        while (this.next[jnode] >= 0) {
            jnode = this.next[jnode];
        }

        this.next[jnode] = inode;
    }
}

function TreeExpand(inode) {
    // alert("in TreeExpand(" + inode + ")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode=" + inode + " is out of range");
        return;
    }

    // expand inode
    this.opened[inode] = 1;

    // update the display
    this.update();
}

function TreeContract(inode) {
    // alert("in TreeContract(" + inode + ")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode=" + inode + " is out of range");
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

function TreeProp(inode, iprop, onoff) {
    // alert("in TreeProp(" + inode + "," + iprop + "," + onoff + ")");

    // validate inputs
    if (inode < 0 || inode >= this.opened.length) {
        alert("inode=" + inode + " is out of range");
        return;
    } else if (onoff != "on" && onoff != "off") {
        alert("onoff=" + onoff + " is not 'on' or 'off'");
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
            thisNode.title = "Toggle Ori off";
        } else {
            thisNode.setAttribute("class", "fakelinkoff");
            thisNode.title = "Toggle Ori on";
        }
    }

    // set property for inode's children
    for (var jnode = inode+1; jnode < this.parent.length; jnode++) {
        if (this.parent[jnode] == inode) {
            this.prop(jnode, iprop, onoff);
        }
    }

    wv.sceneUpd = 1;
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
        newTR.appendChild(newTDcol1);

        var newTexta = doc.createTextNode("");
        newTDcol1.appendChild(newTexta);

        // table data "node"+inode+"col2"
        var newTDcol2 = doc.createElement("TD");
        newTDcol2.setAttribute("id", "node"+inode+"col2");
        if (this.gprim[inode] != "") {
             newTDcol2.className = "fakelinkcmenu";
        }
        //if (this.click[inode] != null) {
        //    newTDcol2.className = "fakelinkcmenu";
        //    if (this.tooltip[inode].length > 0) {
        //        newTDcol2.title = this.tooltip[inode];
        //    }
        //}
        newTR.appendChild(newTDcol2);

        var newTextb = doc.createTextNode(this.name[inode]);
        newTDcol2.appendChild(newTextb);

        var name = this.name[inode].replace(/\u00a0/g, "");

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
}

function TreeUpdate() {
    // alert("in TreeUpdate()");

    var doc = this.document;

    // traverse the nodes using depth-first search
    for (var inode = 1; inode < this.opened.length; inode++) {
        var element = doc.getElementById("node" + inode);

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

        //if (this.click[inode] !== null) {
        //    var myElem = doc.getElementById("node"+inode+"col2");
        //    myElem.onclick = this.click[inode];
        //}

        // set the class of the properties
        if (this.nprop[inode] >= 1) {
            var myElem = doc.getElementById("node"+inode+"col3");
            myElem.onclick = this.cbck1[inode];

            if (this.prop1[inode] == "Viz") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Viz on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Viz off";
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
                        myElem.title = "Toggle Grd on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Grd off";
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
                        myElem.title = "Toggle Trn on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Trn off";
                    }
                }
            } else if (this.prop3[inode] == "Ori") {
                if (this.gprim[inode] != "") {
                    if ((wv.sceneGraph[this.gprim[inode]].attrs & wv.plotAttrs.ORIENTATION) == 0) {
                        myElem.setAttribute("class", "fakelinkoff");
                        myElem.title = "Toggle Ori on";
                    } else {
                        myElem.setAttribute("class", "fakelinkon");
                        myElem.title = "Toggle Ori off";
                    }
                }
            }
        }
    }
}

function TreeClear() {
    // alert("in TreeClear()");
  
    // remove all but the first Node
    this.name.splice(  1);
    this.tooltip.splice(1);
    this.gprim.splice( 1);
    this.parent.splice(1);
    this.child.splice( 1);
    this.next.splice(  1);
    this.nprop.splice( 1);
    this.opened.splice(1);
  
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
}

function Tree(doc, treeId) {
    // alert("in Tree(doc="+doc+", treeId="+treeId+")");

    // remember the document
    this.document = doc;
    this.treeId   = treeId;

    // arrays to hold the Nodes
    this.name    = new Array();
    this.tooltip = new Array();
    this.gprim   = new Array();
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

    // initialize node=0 (the root)
    this.name[   0] = "**root**"
    this.tooltip[0] = "";
    this.gprim[  0] = "*";
    this.parent[ 0] = -1;
    this.child[  0] = -1;
    this.next[   0] = -1;
    this.prop1[ 0]  = "";
    this.cbck1[ 0]  = null;
    this.prop2[ 0]  = "";
    this.cbck2[ 0]  = null;
    this.prop3[ 0]  = "";
    this.cbck3[ 0]  = null;
    this.opened[ 0] = +1;

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
// callback to toggle Viz property
//
function toggleViz(e) {
    // alert("in toggleViz(e="+e+")");

    //if (wv.curMode != 0) {
    //    alert("Command disabled.  Press 'Cancel' or 'OK' first");
    //    return;
    //}

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Viz property
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
            myTree.prop(inode, 1, "on");
        } else {
            myTree.prop(inode, 1, "off");
        }

    //  toggle the Viz property (on all Faces/Edges/Nodes in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col3");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 1, "on");
        } else {
            myTree.prop(inode, 1, "off");
        }
    }

    myTree.update();
}

//
// callback to toggle Trn property
//
function toggleTrn(e) {
    // alert("in toggleTrn(e="+e+")");

    //if (wv.curMode != 0) {
    //    alert("Command disabled.  Press 'Cancel' or 'OK' first");
    //    return;
    //}

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Trn property (on a Face)
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
            myTree.prop(inode, 3, "on");
        } else {
            myTree.prop(inode, 3, "off");
        }

    // toggle the Trn property (on all Faces in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col5");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 3, "on");
            myElem.setAttribute("class", "fakelinkon");
            myElem.title = "Toggle Trn off";
        } else {
            myTree.prop(inode, 3, "off");
            myElem.setAttribute("class", "fakelinkoff");
            myElem.title = "Toggle Trn on";
        }
    }

    myTree.update();
}

//
// callback to toggle Grd property
//
function toggleGrd(e) {
    // alert("in toggleGrd(e="+e+")");

    //if (wv.curMode != 0) {
    //    alert("Command disabled.  Press 'Cancel' or 'OK' first");
    //    return;
    //}

    // get the Tree Node
    var inode = e["target"].id.substring(4);
    inode     = inode.substring(0,inode.length-4);
    inode     = Number(inode);

    // toggle the Grd property
    if (myTree.gprim[inode] != "") {
        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
            myTree.prop(inode, 2, "on");
        } else {
            myTree.prop(inode, 2, "off");
        }

    // toggle the Grd property (on all Faces/Edges/Nodes in this Body)
    } else {
        var myElem = myTree.document.getElementById("node"+inode+"col4");
        if (myElem.getAttribute("class") == "fakelinkoff") {
            myTree.prop(inode, 2, "on");
        } else {
            myTree.prop(inode, 2, "off");
        }
    }

    myTree.update();
};

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


function getMouseDbl(e)
{
    if(document.selection && document.selection.empty) {
        document.selection.empty();
    } else if(window.getSelection) {
        var sel = window.getSelection();
        sel.removeAllRanges();
    }

    wv.centerV = 1;
    postMessage("View being centered");
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


//---------------------------------------------------------------------------//
//
// call-back functions
//
//---------------------------------------------------------------------------//

function wvInitUI()
{

    // set up extra storage for matrix-matrix multiplies
    wv.uiMatrix = new J3DIMatrix4();
    wv.saveMatrix = new J3DIMatrix4(wv.mvMatrix);   // matrix used in save operations

                                   // ui cursor variables
    wv.cursorX   = -1;             // current cursor position
    wv.cursorY   = -1;
    wv.keyPress  = -1;             // last key pressed
    wv.startX    = -1;             // start of dragging position
    wv.startY    = -1;
    wv.button    = -1;             // button pressed
    wv.modifier  =  0;             // modifier (shift,alt,cntl) bitflag
    wv.flying    =  1;             // flying multiplier (do not set to 0)
    wv.offTop    =  0;             // offset to upper-left corner of the canvas
    wv.offLeft   =  0;
    wv.dragging  = false;          // true during drag operation
    wv.picking   =  0;             // keycode of command that turns picking on
    wv.locating  =  0;             // keycode of command that turned locating on
    wv.focus     = [0, 0, 0, 1];   // focus data needed in locating
    wv.debugUI   =  1;             // set to 1 for console messages
    
    wv.gridToggle = true;                // Keep track if all grid states true = on , false = off
    
    wv.transparenceToggle = true;  // Keep track of transparence states
    
    wv.lowerLimit   = -1;          // lower limit in key
    wv.upperLimit   = +1;          // upper limit in key
    
    wv.keyRange = [-1, 1];         // Initial limit range for the color in the key
    
//  wv.pick                        // set to 1 to turn picking on
//  wv.picked                      // sceneGraph object that was picked
//  wv.locate                      // set to 1 to turn on locating
//  wv.sceneGraph                  // pointer to sceneGraph
//  wv.located                     // coordinates (local) that were pointed at
//  wv.centerV                     // set to 1 to put wv into view centering mode
//  wv.sceneUpd                    // set to 1 when scene needs rendering
//  wv.sgUpdate                    // is 1 if the sceneGraph has been updated
//  wv.socketUt.send(text)         // function to send text to server
//  wv.plotAttrs                   // plot attributes   
 
    document.addEventListener('keypress',   getKeyPress,   false);

    var canvas = document.getElementById(wv.canvasID);
    canvas.addEventListener('mousemove',  getCursorXY,     false);
    canvas.addEventListener('mousedown',  getMouseDown,    false);
    canvas.addEventListener('mouseup',    getMouseUp,      false);
    canvas.addEventListener('wheel',      getMouseRoll,    false);
    canvas.addEventListener('dblclick',   getMouseDbl,     false);
    canvas.addEventListener('mouseout',   mouseLeftCanvas, false);
    
    // Add listener to key
    var keycan = document.getElementById(wv.canvasKY);
    keycan.addEventListener('mouseup',    setKeyLimits,    false);
    
}


function wvUpdateUI()
{

    // special code for delayed-picking mode
    if (wv.picking > 0) {

        // if something is picked, post a message
        if (wv.picked !== undefined) {

            // second part of 'g' operation
            if (wv.picking == 103) {
                postMessage("Toggling grid of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.LINES) == 0) {
                            myTree.prop(inode, 2, "on");
                        } else {
                            myTree.prop(inode, 2, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of 't' operation
            } else if (wv.picking == 116) {
                postMessage("Toggling transparency of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.TRANSPARENT) == 0) {
                            myTree.prop(inode, 3, "on");
                        } else {
                            myTree.prop(inode, 3, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of 'v' operation
            } else if (wv.picking == 118) {
                postMessage("Toggling visibility of "+wv.picked.gprim);

                for (var inode = myTree.gprim.length-1; inode >= 0; inode--) {
                    if (myTree.gprim[inode] == wv.picked.gprim) {
                        if ((wv.sceneGraph[myTree.gprim[inode]].attrs & wv.plotAttrs.ON) == 0) {
                            myTree.prop(inode, 1, "on");
                        } else {
                            myTree.prop(inode, 1, "off");
                        }
                        myTree.update();
                        break;
                    }
                }

            // second part of '^' operation
            } else if (wv.picking == 94) {
                var mesg = "Picked: "+wv.picked.gprim;

            //    try {
            //        var attrs = wv.sgData[wv.picked.gprim];
            //        for (var i = 0; i < attrs.length; i+=2) {
            //            mesg = mesg + "\n        "+attrs[i]+"= "+attrs[i+1];
            //        }
            //    } catch (x) {
            //    }
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
    
    // if the tree has not been created but the scene graph (possibly) exists...
    if (wv.sgUpdate == 1 && (wv.sceneGraph !== undefined)) {
    
        // clear previous Nodes from the Tree
        myTree.clear();

        myTree.addNode(0, "Display", "",
                       "Viz", toggleViz,
                       "Grd", toggleGrd);

        // open the Display (by default)
        myTree.opened[1] = 1;

        // ...count number of primitives in the scene graph
        var count = 0;
        for (var gprim in wv.sceneGraph) {
        
            // parse the name
            var matches = gprim.split("|");
        
            // processing when Body is explicitly named: "Bodyname"
            var bodyName = matches[0];

            // process server messages embedded in the name
            if (matches.length >= 2) {
                wvServerMessage(matches.slice(1).join());
            }

            //var ibody = Number(matches[1]);
            //if        (matches[2] == "Face") {
            //    var iface = matches[3];
            //} else if (matches[2] == "Loop") {
            //    var iloop = matches[3];
            //} else if (matches[2] == "Edge") {
            //    var iedge = matches[3];
            //} else {
            //   alert("unknown type: " + matches[2]);
            //    continue;
            //}
        
            // determine if Body does not exists
            var knode = -1;
            for (var jnode = 1; jnode < myTree.name.length; jnode++) {
               if (myTree.name[jnode] == "\u00a0\u00a0"+bodyName) {
                   knode = jnode;
               }
            }
        
            // if Body does not exist, create it and its Face, Loop, and Edge
            //   subnodes now
            var kface, kloop, kedge;
            if (knode < 0) {
               //postMessage("Processing Body " + ibody);
        
               //myTree.addNode(0, "Body " + ibody, "*");
               myTree.addNode(1, "\u00a0\u00a0"+bodyName, gprim,
                              "Viz", toggleViz,
                              "Grd", toggleGrd,
                              "Trn", toggleTrn);
               knode = myTree.name.length - 1;
        
               //myTree.addNode(knode, "__Faces", "*");
               //kface = myTree.name.length - 1;
               //
               //myTree.addNode(knode, "__Loops", "*");
               //kloop = myTree.name.length - 1;
               //
               //myTree.addNode(knode, "__Edges", "*");
               //kedge = myTree.name.length - 1;
        
            // otherwise, get pointers to the face-group and loop-group nodes
            } else {
                kface = myTree.child[knode];
                kloop = kface + 1;
                kedge = kloop + 1;
            }
        
            // make the tree node
            //if        (matches[2] == "Face") {
            //   myTree.addNode(kface, "____face " + iface, gprim);
            //} else if (matches[2] == "Loop") {
            //   myTree.addNode(kloop, "____loop " + iloop, gprim);
            //           } else if (matches[2] == "Edge") {
            //   myTree.addNode(kedge, "____edge " + iedge, gprim);
            //}
        
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
    
        var myKeyPress = String.fromCharCode(wv.keyPress);

        // '?' -- help
        if (wv.keyPress ==  63) {
            postMessage("........................... Viewer Cursor options ............................\n" +
                        "ctrl-h <Home> - initial view             ctrl-f         - front view          \n" +
                        "ctrl-l        - leftside view            ctrl-r         - riteside view       \n" +
                        "ctrl-t        - top view                 ctrl-b         - bottom view         \n" +
                        "ctrl-i <PgUp> - zoom in                  ctrl-o <PgDn>  - zoom out            \n" +
                        "<Left>        - rotate or xlate left     <Rite>         - rotate or xlate rite\n" +
                        "<Up>          - rotate or xlate up       <Down>         - rotate or xlate down\n" +
                        ">             - save view                <              - recall view         \n" +
                        "ctrl-> .      - save view to file        ctrl-< ,      - read view from file \n" +
                        "^ 6           - query object at cursor   @ 2            - get coords. @ cursor\n" +
                        "v             - toggle Viz at cursor     g              - toggle Grd at cursor\n" +
                        "t             - toggle Trn at cursor     o              - toggle Ori at cursor\n" +
                        "                                         * 8 <DblClick> - center view @ cursor\n" +
                        "!             - toggle flying mode       ?              - get help            \n" +
                        "\n" +
                        "........................... DataSet Cursor options ...........................\n" +
                        "n             - next scalar              p              - previous scalar     \n" +
                        "l             - set scalar limits        r              - reverse color map   \n" +
                        "..............................................................................");

        // 'g' -- toggle grid at cursor
        } else if (myKeyPress == 'g' && wv.modifier == 0) {
            wv.picking  = 103;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 't' -- toggle transparency at cursor
        } else if (myKeyPress == 't' && wv.modifier == 0) {
            wv.picking  = 116;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // 'v' -- toggle visibility at cursor
        } else if (myKeyPress == 'v' && wv.modifier == 0) {
            wv.picking  = 118;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '^' -- query at cursor
        } else if (myKeyPress == '^' || myKeyPress == '6') {
            wv.picking  = 94;
            wv.pick     = 1;
            wv.sceneUpd = 1;

        // '@' -- locate at cursor
        } else if (myKeyPress == '@' || myKeyPress == '2') {
            wv.locating = 64;
            wv.locate   = 1;

        // '*' -- center view
        } else if (myKeyPress == '*' || myKeyPress == '8') {
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

        // C-'>' or '.' -- save view to file
        } else if ((wv.keyPress == 46  && wv.modifier == 5) || 
                   (myKeyPress  == "." && wv.modifier == 0)   ) {
            var filename = prompt("Enter view filename to save:", "save.view");
            if (filename !== null) {
                postMessage("Saving view to \"" + filename + "\"");
                browserToServer("saveView|"+filename+"|"+wv.scale+"|"+wv.mvMatrix.getAsArray());
            }

        // C-'<' or ',' -- read view from file
        } else if ((wv.keyPress == 44  && wv.modifier == 5) || 
                   (myKeyPress  == "," && wv.modifier == 0)   ) {
            var filename = prompt("Enter view filename to read:", "save.view");
            if (filename !== null) {
                postMessage("Reading view from \"" + filename + "\"");
                browserToServer("readView|"+filename);
            }

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
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-i' - zoom in (same as <PgUp> without shift)
        } else if ((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) {
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
            cmdHome();
            wv.keyPress = -1;
            return;

        // 'ctrl-r' - riteside view
        } else if ((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) {
            cmdRite();
            wv.keyPress = -1;
            return;

        // 'ctrl-l' - leftside view
        } else if ((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) {
            cmdLeft();
            wv.keyPress = -1;
            return;

        // 'ctrl-t' - top view
        } else if ((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) {
            cmdTop();
            wv.keyPress = -1;
            return;

        // 'ctrl-b' - bottom view
        } else if ((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
                   (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) {
            cmdBotm();
            wv.keyPress = -1;
            return;
      
        // 'n' -- next scalar
        } else if (myKeyPress == 'n' && wv.modifier == 0) {
            browserToServer("nextScalar");
            browserToServer("getRange");
            
            //turnOnColor();
            
        // 'p' -- previous scalar
        } else if (myKeyPress == 'p' && wv.modifier == 0) {
            browserToServer("previousScalar");
            browserToServer("getRange");
            
            //turnOnColor();
            
        // 'l' -- Set scalar limits
        } else if (myKeyPress == 'l' && wv.modifier == 0) {
            
            setKeyLimits(null);
            //turnOnColor();
            
//        // 'C' or 'c' -- Set color map 
//        } else if (wv.keyPress == 67 || wv.keyPress == 99) {
//            
//            var cmap = (prompt("Enter desired color map", cmap))
//            if (cmap !=null) {
//            
//                var colorMap  = "colorMap," + cmap.toLowerCase();
//            
//                //alert(colorMap);
//            
//                var levels = (prompt("Enter desired number of contour levels", levels))
//                if (levels !=null) {
//                
//                    var contourLevel  = "contourLevel," + parseInt(levels);
//                
//                    //alert(contourLevel);
//                
//                    browserToServer(contourLevel);
//                    //turnOnColor();
//                    
//                }
//                
//                browserToServer(colorMap);
//                //turnOnColor();
//            }
//        
        // 'r' -- Reverse/invert the color map  
        } else if (myKeyPress == 'r' && wv.modifier == 0) {
            browserToServer("reverseMap");
            
        // NOP
        } else if (wv.keyPress == 0 && wv.modifier == 0) {

        } else {
            postMessage("'" + String.fromCharCode(wv.keyPress)
                        + "' + modifier '" + wv.modifier 
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
                    var dtheta = Math.atan2(yf, xf)-theta1;
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





function wvServerMessage(text)
{
    // remove trailing NULL
    if (text.charCodeAt(text.length-1) == 0) {
        text = text.substring(0, text.length-1);
    }
   
    
    if (text.startsWith("Range")) {
        
        try {
            
            wv.keyRange= text.split(":")[1].split(",");
            
            document.getElementById("WVkey"   ).hidden = false;
            document.getElementById("CAPSlogo").hidden = true;
            wv.drawKey = 1;
            wv.sceneUpd = 1;

        } catch(err) {
            postMessage(err)
            postParseError(text);
        }
    } else if (text.startsWith("Error")) {
        
        var msgAlert = "Something went wrong... please close the browser." + text
        alert(msgAlert)

    // if it starts with "saveView|" do nothing
    } else if (text.substring(0,9) == "saveView|") {

    // if it starts with "readView|" load view matrix and update display
    } else if (text.substring(0,9) == "readView|") {
        var readViewList = text.split("|");
        if (readViewList.length == 2) {
            postMessage("Could not read view file: \"" + readViewList[1] + "\"");
        } else {
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
        }

    } else {
        
         postMessage(" Server Message: " + text);
    }
    
}


function wvServerDown()
{
    postMessage("The server has terminated or network connection has been lost.");

    // turn the background of the message window pink
    var botm = document.getElementById("brframe");
    botm.style.backgroundColor = "#FFAFAF";
}

//
// callback when "homeButton" is pressed (called by capsViewer.html)
//
var cmdHome = function() {
    wv.mvMatrix.makeIdentity();
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
};

//
// callback when "leftButton" is pressed (called by capsViewer.html)
//
var cmdLeft = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(+90, 0,1,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
};

//
// callback when "riteButton" is pressed (called by capsViewer.html)
//
var cmdRite = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(-90, 0,1,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
};

//
// callback when "botmButton" is pressed (called by capsViewer.html)
//
var cmdBotm = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(-90, 1,0,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
};

//
// callback when "topButton" is pressed (called by capsViewer.html)
//
var cmdTop = function() {
    wv.mvMatrix.makeIdentity();
    wv.mvMatrix.rotate(+90, 1,0,0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale    = 1;
    wv.sceneUpd = 1;
};

//
// callback when "inButton" is pressed (called by capsViewer.html)
//
var cmdIn = function() {
    wv.mvMatrix.load(wv.uiMatrix);
    wv.mvMatrix.scale(2.0, 2.0, 2.0);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale   *= 2.0;
    wv.sceneUpd = 1;
};

//
// callback when "outButton" is pressed (called by capsViewer.html)
//
var cmdOut = function() {
    wv.mvMatrix.load(wv.uiMatrix);
    wv.mvMatrix.scale(0.5, 0.5, 0.5);
    wv.uiMatrix.load(wv.mvMatrix);
    wv.mvMatrix.makeIdentity();
    wv.scale   *= 0.5;
    wv.sceneUpd = 1;
};

// Callback used to put axes on the canvas (called by capsViewer.html)
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
 var alen =  0.25;     // length of axes

 // set up coordinates for axes
 mVal = mv.getAsArray();

 var vertices = new Float32Array(66);
 vertices[ 0] = x;
 vertices[ 1] = y;
 vertices[ 2] = z;
 vertices[ 3] = x + alen*(    mVal[ 0]             );
 vertices[ 4] = y + alen*(    mVal[ 1]             );
 vertices[ 5] = z + alen*(    mVal[ 2]             );
 vertices[ 6] = x + alen*(1.1*mVal[ 0]+0.1*mVal[ 4]);
 vertices[ 7] = y + alen*(1.1*mVal[ 1]+0.1*mVal[ 5]);
 vertices[ 8] = z + alen*(1.1*mVal[ 2]+0.1*mVal[ 6]);
 vertices[ 9] = x + alen*(1.3*mVal[ 0]-0.1*mVal[ 4]);
 vertices[10] = y + alen*(1.3*mVal[ 1]-0.1*mVal[ 5]);
 vertices[11] = z + alen*(1.3*mVal[ 2]-0.1*mVal[ 6]);
 vertices[12] = x + alen*(1.1*mVal[ 0]-0.1*mVal[ 4]);
 vertices[13] = y + alen*(1.1*mVal[ 1]-0.1*mVal[ 5]);
 vertices[14] = z + alen*(1.1*mVal[ 2]-0.1*mVal[ 6]);
 vertices[15] = x + alen*(1.3*mVal[ 0]+0.1*mVal[ 4]);
 vertices[16] = y + alen*(1.3*mVal[ 1]+0.1*mVal[ 5]);
 vertices[17] = z + alen*(1.3*mVal[ 2]+0.1*mVal[ 6]);

 vertices[18] = x;
 vertices[19] = y;
 vertices[20] = z;
 vertices[21] = x + alen*(    mVal[ 4]             );
 vertices[22] = y + alen*(    mVal[ 5]             );
 vertices[23] = z + alen*(    mVal[ 6]             );
 vertices[24] = x + alen*(1.1*mVal[ 4]+0.1*mVal[ 8]);
 vertices[25] = y + alen*(1.1*mVal[ 5]+0.1*mVal[ 9]);
 vertices[26] = z + alen*(1.1*mVal[ 6]+0.1*mVal[10]);
 vertices[27] = x + alen*(1.2*mVal[ 4]             );
 vertices[28] = y + alen*(1.2*mVal[ 5]             );
 vertices[29] = z + alen*(1.2*mVal[ 6]             );
 vertices[30] = x + alen*(1.3*mVal[ 4]+0.1*mVal[ 8]);
 vertices[31] = y + alen*(1.3*mVal[ 5]+0.1*mVal[ 9]);
 vertices[32] = z + alen*(1.3*mVal[ 6]+0.1*mVal[10]);
 vertices[33] = x + alen*(1.2*mVal[ 4]             );
 vertices[34] = y + alen*(1.2*mVal[ 5]             );
 vertices[35] = z + alen*(1.2*mVal[ 6]             );
 vertices[36] = x + alen*(1.2*mVal[ 4]             );
 vertices[37] = y + alen*(1.2*mVal[ 5]             );
 vertices[38] = z + alen*(1.2*mVal[ 6]             );
 vertices[39] = x + alen*(1.2*mVal[ 4]-0.1*mVal[ 8]);
 vertices[40] = y + alen*(1.2*mVal[ 5]-0.1*mVal[ 9]);
 vertices[41] = z + alen*(1.2*mVal[ 6]-0.1*mVal[10]);


 vertices[42] = x;
 vertices[43] = y;
 vertices[44] = z;
 vertices[45] = x + alen*(    mVal[ 8]             );
 vertices[46] = y + alen*(    mVal[ 9]             );
 vertices[47] = z + alen*(    mVal[10]             );
 vertices[48] = x + alen*(1.1*mVal[ 8]+0.1*mVal[ 0]);
 vertices[49] = y + alen*(1.1*mVal[ 9]+0.1*mVal[ 1]);
 vertices[50] = z + alen*(1.1*mVal[10]+0.1*mVal[ 2]);
 vertices[51] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
 vertices[52] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
 vertices[53] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
 vertices[54] = x + alen*(1.3*mVal[ 8]+0.1*mVal[ 0]);
 vertices[55] = y + alen*(1.3*mVal[ 9]+0.1*mVal[ 1]);
 vertices[56] = z + alen*(1.3*mVal[10]+0.1*mVal[ 2]);
 vertices[57] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
 vertices[58] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
 vertices[59] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
 vertices[60] = x + alen*(1.1*mVal[ 8]-0.1*mVal[ 0]);
 vertices[61] = y + alen*(1.1*mVal[ 9]-0.1*mVal[ 1]);
 vertices[62] = z + alen*(1.1*mVal[10]-0.1*mVal[ 2]);
 vertices[63] = x + alen*(1.3*mVal[ 8]-0.1*mVal[ 0]);
 vertices[64] = y + alen*(1.3*mVal[ 9]-0.1*mVal[ 1]);
 vertices[65] = z + alen*(1.3*mVal[10]-0.1*mVal[ 2]);

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
 gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
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
}


//function turnOnColor() {
//    
//    for (var gprim in wv.sceneGraph) { // Loop through scene primatives
            
        //wv.sceneGraph[gprim].attrs |=  wv.plotAttrs.SHADING;
        //wv.sceneGraph[gprim].attrs |=  wv.plotAttrs.LINES;
//    }    
//}

// Callback when the mouse is pressed in key window
function setKeyLimits(e) {
    
    try {
        wv.lowerLimit = wv.key.bottom;
        wv.upperLimit = wv.key.top;
    } catch(err) {
        postMessage("Error: no key set yet!")
        return 
    } 
    
    var msg = "Enter upper limit (max = " + wv.keyRange[1] + ")" 
    var tempup = prompt(msg, wv.upperLimit);
    if (tempup == null) {
        return; 
    } 
    
    if (isNaN(tempup)) {
        alert("Upper limit must be a number");
        return;
    }
    
    // Get new limits
    msg = "Enter lower limit (min = " + wv.keyRange[0] + ")" 
    var templo = prompt(msg, wv.lowerLimit);
    
    if (templo == null) {
        return; 
    } 
    
    if (isNaN(templo)) {
     alert("Lower limit must be a number");
     return;
    }


    if (Number(tempup) < Number(templo)) {
        alert("Upper limit must be greater than lower limit");
        return;
    }

    if (templo != wv.lowerLimit || tempup != wv.uplimit) {
        wv.lowerLimit = templo;
        wv.upperLimit = tempup;

        // send the limits back to the server
        var limits = "limits,"+ wv.lowerLimit.toString() + "," + wv.upperLimit.toString();
        browserToServer(limits);
    }
    
    wv.drawKey = 1;
}
