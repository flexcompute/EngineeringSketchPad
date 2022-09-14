const svg       = d3.select("#viz-area");
const svgWidth  = svg.style("width");
const svgHeight = svg.style("height");

const colors = [
    "#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"
];

const strokes = [
    "5,5", "15,7,5,3,5,7", "10,10"
];

const data = JSON.parse(JSON.stringify(dataJSON));

// process data
for (let i=0; i < data.aims.length; i++) {
    const aim = data.aims[i];
    aim["valIn"] = [];
    aim["valOut"] = [];
    aim["geomIn"] = [];
    aim["geomOut"] = [];
    aim["coords"] = {
        "x1": 0,
        "x2": 0,
        "y1": 0,
        "y2": 0
    };
    const colorIndex = i % colors.length;
    const strokeIndex = i % strokes.length;
    aim["outColor"] = colors[colorIndex];
    aim["strokeArray"] = strokes[strokeIndex];
}

// push link in and out data to each node
for (let i=0; i < data.valLinks.length; i++) {
    const link = data.valLinks[i];
    const source = data.aims.find(x => x.id === link.source);
    const target = data.aims.find(x => x.id === link.target);
    if (source.valOut.includes(link.target)) {
        const used = data.valLinks.find(el => (el.source === link.source) && (el.target === link.target) && el.include);
        used.data.push({"sourceVar": link.data[0].sourceVar, "targetVar": link.data[0].targetVar});
    } else {
        source.valOut.push(link.target);
        target.valIn.push(link.source);
        link["include"] = true;
    }
}

for (let i=0; i < data.geomLinks.length; i++) {
    const link = data.geomLinks[i];
    const source = data.aims.find(x => x.id === link.source);
    const target = data.aims.find(x => x.id === link.target);
    if (source.geomOut.includes(link.target)) {
        const used = data.geomLinks.find(el => (el.source === link.source) && (el.target === link.target) && el.include);
        used.data.push({"bound": link.data[0].bound, "sourceVar": link.data[0].sourceVar, "targetVar": link.data[0].targetVar})
    } else {
        source.geomOut.push(link.target);
        target.geomIn.push(link.source);
        link["include"] = true;
    }
}

const ranksep      = 100;
const svgBoxHeight = 60;
const svgBoxWidth  = 180;

const graph = new dagre.graphlib.Graph();
graph.setGraph({
    'marginx': 50,
    'marginy': 50,
    'acyclicer': 'greedy',
    'ranker': 'longest-path', //network-simplex, tight-tree
    'ranksep': ranksep
});
graph.setDefaultEdgeLabel(() => { return {}; });

for (let i=0; i < data.aims.length; i++) {
    const aim = data.aims[i];
    // const width = aim.id.length * 15
    graph.setNode(aim["id"], { label: aim.id, width: svgBoxWidth, height: svgBoxHeight });
}

for (let i=0; i < data.valLinks.length; i++) {
    const edge = data.valLinks[i];
    if (edge["include"]) {
        graph.setEdge(edge.source, edge.target);
    }
}

dagre.layout(graph);

const graphWidth = graph.graph().width;
const graphHeight = graph.graph().height;
const viewBox = '0 0 ' + graphWidth + ' ' + graphHeight;
svg.attr('viewBox', viewBox);

const textDiv = document.getElementById('textarea');

// scaling function for d3 between percent and pixel
const x = d3.scaleLinear()
      .domain([0, 100])         
      .range([0, svgWidth]);       

const y = d3.scaleLinear()
      .domain([0, 100])
      .range([0, svgHeight]);

const aimDiv   = d3.select('div#aim-hover');
const aimDivJs = document.getElementById('aim-hover');

// Draw a box for each AIM
graph.nodes().forEach((v) => {
    const node = graph.node(v);
    const boxWidth = node.width;
    const boxHeight = node.height;
    const box_x = node.x - boxWidth/2;
    const box_y = node.y - boxHeight/2;
    const data_aim = data.aims.find(x => x.id == v);
    const inVars = data_aim.inVars;
    const outVars = data_aim.outVars;
    const dynOutVars = data_aim.dynOutVars;
    data_aim.coords.x1 = box_x;
    data_aim.coords.x2 = box_x + boxWidth;
    data_aim.coords.y1 = box_y;
    data_aim.coords.y2 = box_y + boxHeight;
    const boxColor = data_aim.dirty ? "#FF0000" : "#111111";
    const aimGroup = svg.append("g");
    const outColor = (data_aim.valOut.length > 0 || data_aim.geomOut.length > 0) ? data_aim.outColor : boxColor;
    aimGroup.append("rect")
        .attr("x", box_x)
        .attr("y", box_y)
        .attr("width", boxWidth)
        .attr("height", boxHeight)
        .attr("fill", boxColor)
        .attr("stroke", outColor)
        .attr("stroke-width", 2);
    aimGroup.append("text")
        .attr("x", box_x + boxWidth/2)
        .attr("y", box_y + boxHeight/2)
        .attr("stroke", "white")
        .attr("fill", "white")
        .style("font-family", '"Open Sans", sans-serif')
        .style("text-anchor", "middle")
        .style("alignment-baseline", "middle")
        .style("letter-spacing", "0.1em")
        .text(node.label);
    aimGroup.append("rect") // invisible box for mouseover event
        .attr("x", box_x)
        .attr("y", box_y)
        .attr("width", boxWidth)
        .attr("height", boxHeight)
        .attr("opacity", "0");
    aimGroup.on('mouseover', function (event) {
        // can't use arrow function without changing 'this' to global
        d3.select(this).transition()
            .duration('50')
            .attr('opacity', '.85');
        aimDiv.selectAll('tr.line').remove();
        aimDiv.transition()
            .duration('50')
            .style('opacity', 1);
        aimDiv
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 15) + 'px');
        const hoverTitle = data_aim.dirty ? node.label + " (dirty)" : node.label;
        aimDiv.select('#hover-title').text(hoverTitle);

        const inTable = aimDiv.select('table#aim-table-in');
        const outTable = aimDiv.select('table#aim-table-out');
        for (let i=0; i<inVars.length; i++) {
            const varName = inVars[i].name;
            // take tuple value and stringify if type object
            // stringify options: round if number to four digits, indent using 2 spaces
            let varVal;
            if (typeof inVars[i].value === 'object' && inVars[i].value !== null) {
                varVal = JSON.stringify(inVars[i].value, 
                                        (key, val) => {
                                            return val.toFixed ? Number(val.toFixed(4)) : val
                                        }, 
                                        2)
            } else {
                varVal = inVars[i].value;
            }
            const tableRow = inTable.append('tr').attr('class', 'line');
            tableRow.append('td').text(varName)
            const valueCell = tableRow.append('td');
            valueCell.append('pre').text(varVal);
        }
        for (let i=0; i<outVars.length; i++) {
            const varName = outVars[i].name;
            let varVal;
            if (typeof outVars[i].value === 'object' && outVars[i].value !== null) {
                varVal = JSON.stringify(outVars[i].value, 
                                        (key, val) => {
                                            return val.toFixed ? Number(val.toFixed(4)) : val
                                        }, 
                                        2)
            } else {
                varVal = outVars[i].value;
            }
            const tableRow = outTable.append('tr').attr('class', 'line');
            tableRow.append('td').text(varName);
            const valueCell = tableRow.append('td');
            valueCell.append('pre').text(varVal);
        }
        if (dynOutVars.length > 0) {
            const dividerRow = outTable.append('tr').attr('class', 'line').attr('id', 'divider');
            dividerRow.append('td').text("Dynamic Outputs");
            dividerRow.append('td');
        }
        for (let i=0; i<dynOutVars.length; i++) {
            const varName = dynOutVars[i].name;
            let varVal;
            if (typeof dynOutVars[i].value === 'object' && dynOutVars[i].value !== null) {
                varVal = JSON.stringify(dynOutVars[i].value, 
                                        (key, val) => {
                                            return val.toFixed ? Number(val.toFixed(4)) : val
                                        }, 
                                        2)
            } else {
                varVal = dynOutVars[i].value;
            }
            const tableRow = outTable.append('tr').attr('class', 'line');
            tableRow.append('td').text(varName);
            const valueCell = tableRow.append('td');
            valueCell.append('pre').text(varVal);
        }
    })
        .on('mouseout', function (event) {
            d3.select(this).transition()
                .duration('50')
                .attr('opacity', '1');
            if (!aimDivJs.matches(":hover")) {
                aimDiv.transition()
                    .duration('50')
                    .style('opacity', 0)
                    .on('end', () => {
                        aimDiv.selectAll('tr.line').remove();
                    });
            }
            
        })
        .on('click', function(event) {
            const header = data_aim.id + " Variables";
            
            const textBlock = document.createElement("div")
            textBlock.className = "text-block"
            
            const heading = document.createElement("h3");
            const headingText = document.createTextNode(header);
            heading.appendChild(headingText);
            textBlock.appendChild(heading);
            
            const template = document.getElementById("textbox-table-template");
            const clone = template.content.cloneNode(true);
            
            const inTable = clone.querySelector("#textbox-in-table tbody");
            const outTable = clone.querySelector("#textbox-out-table tbody");
            
            for (let i=0; i<inVars.length; i++) {
                const varName = inVars[i].name;
                let varVal;
                if (typeof inVars[i].value === 'object' && inVars[i].value !== null) {
                    varVal = JSON.stringify(inVars[i].value, 
                                            (key, val) => {
                                                return val.toFixed ? Number(val.toFixed(4)) : val
                                            }, 
                                            2)
                } else {
                    varVal = inVars[i].value;
                }
                const tr = document.createElement("tr")
                const td1 = document.createElement('td')
                const td2 = document.createElement('td')
                const pre = document.createElement('pre')
                td1.textContent = varName;
                pre.textContent = varVal;
                td2.appendChild(pre);
                tr.appendChild(td1);
                tr.appendChild(td2);
                inTable.appendChild(tr);
            }
            for (let i=0; i<outVars.length; i++) {
                const varName = outVars[i].name;
                let varVal;
                if (typeof outVars[i].value === 'object' && outVars[i].value !== null) {
                    varVal = JSON.stringify(outVars[i].value, 
                                            (key, val) => {
                                                return val.toFixed ? Number(val.toFixed(4)) : val
                                            }, 
                                            2)
                } else {
                    varVal = outVars[i].value;
                }
                let varDeriv = null;
                if (outVars[i].deriv !== null) {
                    varDeriv = JSON.stringify(outVars[i].deriv,
                                              (key, val) => {
                                                  return val.toFixed ? Number(val.toFixed(4)) : val
                                              }, 
                                              2);
                }
                const tr = document.createElement("tr")
                const td1 = document.createElement('td')
                const td2 = document.createElement('td')
                const td3 = document.createElement('td')
                const pre2 = document.createElement('pre')
                const pre3 = document.createElement('pre');
                td1.textContent = varName;
                pre2.textContent = varVal;
                if (varDeriv) {
                    pre3.textContent = varDeriv;
                }
                td2.appendChild(pre2);
                td3.appendChild(pre3);
                tr.appendChild(td1);
                tr.appendChild(td2);
                tr.appendChild(td3);
                outTable.appendChild(tr);
            }
            if (dynOutVars.length > 0) {
                // const dividerRow = outTable.append('tr').attr('class', 'line').attr('id', 'divider');
                // dividerRow.append('td').text("Dynamic Variables");
                const tr = document.createElement("tr");
                tr.id = "divider";
                const td1 = document.createElement("td");
                const td2 = document.createElement("td");
                const td3 = document.createElement("td");
                td1.textContent = "Dynamic Outputs";
                tr.appendChild(td1);
                tr.appendChild(td2);
                tr.appendChild(td3);
                outTable.appendChild(tr);
            }
            for (let i=0; i<dynOutVars.length; i++) {
                const varName = dynOutVars[i].name;
                let varVal;
                if (typeof dynOutVars[i].value === 'object' && dynOutVars[i].value !== null) {
                    varVal = JSON.stringify(dynOutVars[i].value, 
                                            (key, val) => {
                                                return val.toFixed ? Number(val.toFixed(4)) : val
                                            }, 
                                            2)
                } else {
                    varVal = dynOutVars[i].value;
                }
                let varDeriv = null;
                if (dynOutVars[i].deriv !== null) {
                    varDeriv = JSON.stringify(dynOutVars[i].deriv,
                                              (key, val) => {
                                                  return val.toFixed ? Number(val.toFixed(4)) : val
                                              }, 
                                              2);
                }
                const tr = document.createElement("tr");
                const td1 = document.createElement('td');
                const td2 = document.createElement('td');
                const td3 = document.createElement('td');
                const pre2 = document.createElement('pre');
                const pre3 = document.createElement('pre');
                td1.textContent = varName;
                pre2.textContent = varVal;
                if (varDeriv) {
                    pre3.textContent = varDeriv;
                }
                td2.appendChild(pre2);
                td3.appendChild(pre3);
                tr.appendChild(td1);
                tr.appendChild(td2);
                tr.appendChild(td3);
                outTable.appendChild(tr);
            }
            textBlock.appendChild(clone)
            textDiv.appendChild(textBlock)  
            
            textBlock.scrollIntoView();
        });
});

const markerBoxWidth = 10;
const markerBoxHeight = 10;
const refX = markerBoxWidth;
const refY = markerBoxHeight / 2;
const markerWidth = markerBoxWidth / 2;
const markerHeight = markerBoxHeight / 2;
const arrowPoints = [[0, 0], [0, markerBoxHeight], [markerBoxWidth, markerBoxHeight/2], [0, 0]];

svg
    .append('defs')
    .append('marker')
    .attr('id', 'arrow')
    .attr('viewBox', [0, 0, markerBoxWidth, markerBoxHeight])
    .attr('refX', refX)
    .attr('refY', refY)
    .attr('markerWidth', markerBoxWidth)
    .attr('markerHeight', markerBoxHeight)
    .attr('orient', 'auto')
    .append('path')
    .attr('d', d3.line()(arrowPoints))
    .attr('stroke', 'black');

const linkDiv = d3.select('div#link-hover');
const linkDivJs = document.getElementById('link-hover');

graph.edges().forEach((e, i) => {
    const edge = graph.edge(e);
    const points = edge.points;
    const source = data.aims.find(x => x.id == e.v);
    const target = data.aims.find(x => x.id == e.w);
    const dataEdge = data.valLinks.find(x => x.source == e.v && x.target == e.w && x.include);
    const path = points.map((v) => [v.x, v.y])
    // create colored arrowhead marker
    const arrowId = 'arrow' + String(i);
    const arrowUrl = 'url(#' + arrowId + ')';
    const linkGroup = svg.append("g");
    linkGroup
        .append('defs')
        .append('marker')
        .attr('id', arrowId)
        .attr('viewBox', [0, 0, markerBoxWidth, markerBoxHeight])
        .attr('refX', refX)
        .attr('refY', refY)
        .attr('markerWidth', markerBoxWidth)
        .attr('markerHeight', markerBoxHeight)
        .attr('orient', 'auto')
        .append('path')
        .attr('d', d3.line()(arrowPoints))
        .attr('stroke', source.outColor)
        .attr('fill', source.outColor)
        .attr('pointer-events', 'all');
    
    linkGroup.append('path')
        .attr('d', d3.line()(path))
        .attr('stroke', source.outColor)
        .attr('stroke-width', 1.5)
        .attr('marker-end', arrowUrl)
        .attr('fill', 'none')
        .attr('class', 'visible-line-path');

    linkGroup.append('path')
        .attr('d', d3.line()(path))
        .attr('stroke', source.outColor)
        .attr('stroke-width', 20)
        .attr('marker-end', arrowUrl)
        .attr('fill', 'none')
        .attr('opacity', '0');
  
  
    linkGroup.on('mouseover', function (event) {
        linkDiv.selectAll('ul').html("");
        d3.select(this).transition()
            .duration('50')
            .attr('opacity', '.85');
        d3.select(this).select('path.visible-line-path').transition()
            .duration('50')
            .attr('stroke-width', 2)
        linkDiv.transition()
            .duration('50')
            .style('opacity', 1);
        linkDiv
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 15) + 'px');
        const linkTitle = source.id + ' to ' + target.id;
        linkDiv.select('#hover-title').text(linkTitle);
        const inVar = linkDiv.select('ul#link-left-ul');
        const outVar = linkDiv.select('ul#link-right-ul');
        for (let j=0; j<dataEdge.data.length; j++) {
            inVar.append('li').text(dataEdge.data[j].sourceVar);
            outVar.append('li').text(dataEdge.data[j].targetVar);
        }
    })
        .on('mouseout', function(d, i) {
            d3.select(this).transition()
                .duration('50')
                .attr('opacity', '1');
            d3.select(this).select('path.visible-line-path').transition()
                .duration('50')
                .attr('stroke-width', 1.5)
            linkDiv.transition()
                .duration('50')
                .style('opacity', 0)
                .on('end', () => {
                    linkDiv.selectAll('ul').html("");
                });
            
        })
        .on('click', function() {
            const header = "Link from " + source.id + " to " + target.id;

            const textBlock = document.createElement("div")
            textBlock.className = "text-block"
  
            const heading = document.createElement("h3");
            const headingText = document.createTextNode(header);
            heading.appendChild(headingText);
            textBlock.appendChild(heading);

            for (let j=0; j < dataEdge.data.length; j++) {
                const src = dataEdge.data[j].sourceVar;
                const tgt = dataEdge.data[j].targetVar;
                const text = "Out: " + src + "\r\n In: " + tgt;
                const p = document.createElement("p");
                p.textContent = text;
                p.className = "textbox-link"
                textBlock.appendChild(p);
            }

            textDiv.appendChild(textBlock)  
            
            textBlock.scrollIntoView();
        });
});

const geomDiv = d3.select('div#geom-hover');

const usedGeomLinks = data.geomLinks.filter(edge => edge.include);

for (let i=0; i < usedGeomLinks.length; i++) {
    const link = usedGeomLinks[i];
    const source = data.aims.find(x => x.id === link.source);
    const target = data.aims.find(x => x.id === link.target);
    const startX = source.coords.x2;
    const startY = 0.5 * (source.coords.y1 + source.coords.y2);
    const endX = target.coords.x1;
    const endY = 0.5 * (target.coords.y1 + target.coords.y2);
    let path;
    if (endX > startX) { // if left of 2nd box is to the right
        midX = 0.5 * startX + 0.5 * endX;
        path = [[startX, startY], [midX, startY], [midX, endY], [endX, endY]];
    } else if (endX < startX) { // if left of 2nd box is to the left
        const x1 = startX + 30;
        const x2 = endX - 30;
        const offsetY = startY + svgBoxHeight/2 + ranksep/3;
        path = [[startX, startY], [x1, startY], [x1, offsetY], [x2, offsetY], [x2, endY], [endX, endY]];
    } 

  // create colored arrowhead marker
    const arrowId = 'arrowG' + String(i);
    const arrowUrl = 'url(#' + arrowId + ')';
    const linkGroup = svg.append("g");
    linkGroup
        .append('defs')
        .append('marker')
        .attr('id', arrowId)
        .attr('viewBox', [0, 0, markerBoxWidth, markerBoxHeight])
        .attr('refX', refX)
        .attr('refY', refY)
        .attr('markerWidth', markerBoxWidth)
        .attr('markerHeight', markerBoxHeight)
        .attr('orient', 'auto')
        .append('path')
        .attr('d', d3.line()(arrowPoints))
        .attr('stroke', source.outColor)
        .attr('fill', source.outColor);
    linkGroup.append('path')
        .attr('d', d3.line()(path))
        .attr('stroke', source.outColor)
        .attr('stroke-width', 1.5)
        .attr('stroke-dasharray', source.strokeArray)
        .attr('marker-end', arrowUrl)
        .attr('fill', 'none')
        .attr('class', 'visible-line-path');
    linkGroup.append('path')
        .attr('d', d3.line()(path))
        .attr('stroke', source.outColor)
        .attr('stroke-width', 10)
        .attr('marker-end', arrowUrl)
        .attr('fill', 'none')
        .attr('opacity', '0');
    
    linkGroup.on('mouseover', function (event) {
    d3.select(this).transition()
            .duration('50')
            .attr('opacity', '.85');
        d3.select(this).select('path.visible-line-path').transition()
            .duration('50')
            .attr('stroke-width', 2)
        geomDiv.transition()
            .duration('50')
            .style('opacity', 1);
        geomDiv
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 15) + 'px');
        const linkTitle = source.id + ' to ' + target.id;
        geomDiv.select('#hover-title').text(linkTitle)
        const bound = geomDiv.select('ul#bound-ul');
        const inVar = geomDiv.select('ul#geom-left-ul');
        const outVar = geomDiv.select('ul#geom-right-ul');
        for (let j=0; j < link.data.length; j++) {
            bound.append('li').text(link.data[j].bound);
            inVar.append('li').text(link.data[j].sourceVar);
            outVar.append('li').text(link.data[j].targetVar);
        }
    })
        .on('mouseout', function(d, i) {
            d3.select(this).transition()
                .duration('50')
                .attr('opacity', '1');
            d3.select(this).select('path.visible-line-path').transition()
                .duration('50')
                .attr('stroke-width', 1.5)
            geomDiv.transition()
                .duration('50')
                .style('opacity', 0)
                .on('end', () => {
                    geomDiv.selectAll('ul').html("");
                });
            
        })
        .on('click', function() {
            const header = "Link from " + source.id + " to " + target.id;
            const textBlock = document.createElement("div")
            textBlock.className = "text-block"
            
            const heading = document.createElement("h3");
            const headingText = document.createTextNode(header);
            heading.appendChild(headingText);
            textBlock.appendChild(heading);
            
            for (let j=0; j < link.data.length; j++) {
                const bound = link.data[j].bound;
                const src = link.data[j].sourceVar;
                const tgt = link.data[j].targetVar;
                const text = "Bound: " + bound + "\r\n Out: " + src + "\r\n In: " + tgt;
                const p = document.createElement("p");
                p.className = "textbox-link"
                p.textContent = text;
                textBlock.appendChild(p);
            }
            
            textDiv.appendChild(textBlock)  
            textBlock.scrollIntoView();
        });
}

const controlButtons = false
const colorButton = document.getElementById('toggle-color');
const strokeButton = document.getElementById('toggle-stroke');

if (controlButtons) {
    const lineColors = svg.selectAll('path').nodes().map(path => path.getAttribute('stroke'));
    const lineMarkers = svg.selectAll('path').nodes().map(path => path.getAttribute('marker-end'));
    const lineStrokes = svg.selectAll('path').nodes().map(path => path.getAttribute('stroke-dasharray'));


    colorButton.onclick = () => {
        const currentColors = svg.selectAll('path').nodes().map(path => path.getAttribute('stroke'));
        const currentMarkers = svg.selectAll('path').nodes().map(path => path.getAttribute('marker-end'));
        svg.selectAll('path')
            .attr('stroke', (d, i) => {
                if (currentColors[i] == 'black') {
                    return lineColors[i];
                } 
                return 'black';
            })
            .attr('marker-end', (d, i) => {
                if (currentMarkers[i] == 'url(#arrow)') {
                    return lineMarkers[i];
                }
                return 'url(#arrow)';
            })
    }


    strokeButton.onclick = () => {
        const currentStrokes = svg.selectAll('path').nodes().map(path => path.getAttribute('stroke-dasharray'));
        svg.selectAll('path')
            .attr('stroke-dasharray', (d, i) => {
                if (currentStrokes[i] == '1,0') {
                    return lineStrokes[i];
                } 
                return '1,0';
            })
    }
} else {
    colorButton.style.display = 'none';
    strokeButton.style.display = 'none';
}


const BORDER_SIZE = 4;
const textAreaDiv = document.getElementById('textarea-container');
const menuDiv = document.getElementById('menu-div');
const svgAreaDiv = document.getElementById('svg-div');
let m_pos;

function resize(e) {
    const dy = m_pos - e.y;
    m_pos = e.y;
    textAreaDiv.style.height = (parseInt(getComputedStyle(textAreaDiv, '').height) + dy) + "px";
    const newSvgHeight = parseInt(textAreaDiv.style.height) + parseInt(getComputedStyle(menuDiv).height);
    const svgHeightStr = String(newSvgHeight);
    svgAreaDiv.style.height = "calc(100vh - " + svgHeightStr + "px)";
}

textAreaDiv.addEventListener("mousedown", function(e){
    if (e.offsetY < BORDER_SIZE) {
        m_pos = e.y;
        document.addEventListener("mousemove", resize, false);
    }
}, false);

document.addEventListener("mouseup", function(){
    document.removeEventListener("mousemove", resize, false);
}, false);

document.addEventListener("mouseup", function(){
    document.removeEventListener("mousemove", resize, false);
}, false);

aimDivJs.addEventListener("mouseleave", () => {
    aimDiv.transition()
        .duration('50')
        .style('opacity', 0)
        .on('end', () => {
            aimDiv.selectAll('tr.line').remove();
        });
});
