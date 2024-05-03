"use strict"

//--------------vis network-------------
var nodes = null;
var edges = null;
var network = null;

var WIDTH_SCALE = 2,
  GREEN = "#1cc88a",
  RED = "#e74a3b",
  GRAY = "#858796",
  BLACK = "#2B1B17",
  YELLOW = "#f6c23e",
  LIGHTBLUE = "#36b9cc",
  BLUE = "#4e73df";

var microgridData = {
    nodes: [
        {id: 1, label: "Controller 1", group: "Node", value: 10},
        {id: 2, label: "Controller 2", group: "Node", value: 8},
        {id: 3, label: "Controller 3", group: "Node", value: 6},
        // Add other nodes here...
    ],
    edges: [
        {from: 1, to: 2, length: 100, width: WIDTH_SCALE * 6, label: "100"},
        {from: 1, to: 3, length: 150, width: WIDTH_SCALE * 4, label: "150"},
        // Add other edges here...
    ]
};

// Called when the Visualization API is loaded.
function draw(microgridData) {
  nodes = new vis.DataSet(microgridData.nodes);
  edges = new vis.DataSet(microgridData.edges);


// -----------------legend--------------------------
  var mynetwork = document.getElementById("mynetwork");
  var x = -mynetwork.clientWidth / 2 + 50;
  var y = -mynetwork.clientHeight / 2 + 50;
  var step = 70;
  // Legend nodes
  /*var legendNodes = [
    {id: 1000, x: x, y: y, label: "Generator", group: "generator", value: 1, fixed: true, physics: false},
    {id: 1001, x: x, y: y + step, label: "Wind Turbine", group: "windTurbine", value: 1, fixed: true, physics: false},
    {id: 1002, x: x, y: y + 2 * step, label: "Solar Panel", group: "solarPanel", value: 1, fixed: true, physics: false},
    {id: 1003, x: x, y: y + 3 * step, label: "Battery", group: "batteryStorage", value: 1, fixed: true, physics: false},
    {id: 1004, x: x, y: y + 4 * step, label: "Load", group: "criticalLoad", value: 1, fixed: true, physics: false},
    {id: 1005, x: x, y: y + 5 * step, label: "Node", group: "Node", value: 1, fixed: true, physics: false}
  ];
  // Add legend nodes to the nodes DataSet
  legendNodes.forEach(node => nodes.add(node));*/

// -----------------create a network-----------------
  var container = document.getElementById("mynetwork");
  var data = {
    nodes: nodes,
    edges: edges,
  };
  var options = {
    autoResize: true,
    nodes: {
      /*scaling: {
        min: 16,
        max: 32,
      },*/
    },
    edges: {
      color: GRAY,
      smooth: false,
      font: {
            align: 'top', // 'horizontal', 'top', 'middle', or 'bottom'
      },
      arrows: {
          to: { enabled: true, scaleFactor: 1, type: 'arrow' },
      }
    },
    physics: {
        enabled: true,
        barnesHut: {
            gravitationalConstant: -35000,
            centralGravity: 1,
            springLength: 1,
            springConstant: 0.01
        }
    },
    interaction: {
        dragNodes: true,
        dragView: true
    },
    groups: {
      generator: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\ue932',
            size: 50,
            color: GREEN,
        }
      },
      windTurbine: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\uec0c',
            size: 50,
            color: LIGHTBLUE,
        }
      },
      solarPanel: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\uec0f',
            size: 50,
            color: YELLOW,
        }
      },
      batteryStorage: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\ue1a3',
            size: 50,
            color: GREEN,
        }
      },
      criticalLoad: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\uea40',
            size: 50,
            color: BLUE,
        }
      },
      Node: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\uef4a',
            size: 50,
            color: BLUE,
        }
      },
    },
  };
  network = new vis.Network(container, data, options);
  // Fit the network once the stabilization is done
  // Function to check and adjust the view
    

    // Listen to dragEnd and zoom events to adjust the view
    //network.on("dragEnd", function(params) {
    //    checkAndAdjustView();
    //});

    //network.on("zoom", function(params) {
    //    checkAndAdjustView();
    //});

    // Add click event listener for nodes
    network.on("click", function (params) {
        var nodeInfoPanel = document.getElementById("nodeInfoPanel");
        var nodeInfoContent = document.getElementById("nodeInfoContent");
        // Hide the panel initially on every click, then show as needed
        nodeInfoPanel.style.display = "none";

        // When click on nodes
        if (params.nodes.length > 0) {
            var nodeId = params.nodes[0]; // Get the first clicked node ID
            var nodeData = nodes.get(nodeId); // Retrieve the node data from the DataSet

            // Construct the content to display
            var content = `<p>ID: ${nodeData.id}</p>
            `;
            if (edgeData.edgetype == "Node"){
                content += `<p>Node Type: Bus</p>`;
            }else{
                content +=`<p>Node Type: ${nodeData.group}</p>`;
            }

            if (nodeData.group =="generator"){

                content += `<p>Generation Power: ${Number(nodeData.value.toFixed(2))}W</p>`;
            }
            else{
                content += `<p>Nominal Voltage: ${Number(nodeData.value.toFixed(2))}V</p>`;
            }
            // Update the info panel with the node data
            nodeInfoContent.innerHTML = content;

            // Show the info panel
            nodeInfoPanel.style.display = "block";
            nodeInfoPanel.style.top = '80px';
        }
        // When an edge was clicked
        else if (params.edges.length > 0) {
            var edgeId = params.edges[0];
            var edgeData = edges.get(edgeId);


            var content = `<p>Edge Type: ${edgeData.edgetype}</p>
                           <p>Length: ${Number(edgeData.length.toFixed(2))}km</p>`;

            if (edgeData.edgetype == "Transformer"){
                content += `<p>Primary Voltage: ${Number(edgeData.primaryVoltage.toFixed(2))}V</p>`;
                content += `<p>Secondary Voltage: ${Number(edgeData.secondaryVoltage.toFixed(2))}V</p>`;
                content += `<p>Power Rating: ${Number(edgeData.powerRating.toFixed(2))}kVA</p>`;
            }
            else if (edgeData.edgetype == "Overhead lines" || edgeData.edgetype == "Underground lines"){
                content += `<p>Power Flow: ${Number(edgeData.powerFlow.toFixed(2))}W</p>`;
            }


            nodeInfoContent.innerHTML = content;
            nodeInfoPanel.style.display = "block";
            nodeInfoPanel.style.top = '80px';
        }

    });


}

//window.addEventListener("load", () => {
//  draw(microgridData);
//});



function updateGraphBasedOnHour() {
    var hour = document.getElementById('hourSlider').value;
    document.getElementById('hourDisplay').innerText = hour;

    // Find max and min value for multi outputs
    let maxValue = 0;
    let minValue = Number.MAX_SAFE_INTEGER;
    edges.forEach((edge) => {
        if (edge.multiOutputs && edge.multiOutputs.length > hour) {
            maxValue = Math.max(maxValue, edge.multiOutputs[hour]);
            minValue = Math.min(minValue, edge.multiOutputs[hour]);
        }
    });

    // Define the maximum and minimum widths
    const maxWidth = 7;
    const minWidth = 0.5;

    // Iterate over all edges and update their width based on the selected hour
    edges.forEach((edge, id) => {
        if (edge.multiOutputs && edge.multiOutputs.length > hour) {
            var output = edge.multiOutputs[hour];
            // Normalize the width
            var abs_output = Math.abs(output);
            var normalizedWidth = minWidth + (abs_output - minValue) / (maxValue - minValue) * (maxWidth - minWidth);
            normalizedWidth = Math.max(minWidth, Math.min(maxWidth, normalizedWidth)); // Ensure within bounds
            //console.log('Normalized Width for edge', id, 'at hour', hour, 'is', normalizedWidth);
            var fromnode = edge.from;
            var tonode = edge.to;
            var direction = edge.direction;
            if(output < 0 && edge.direction == 'positive'){
                fromnode = edge.to;
                tonode = edge.from;
                direction = "negative";
                console.log('Reverse Direction for Edge', id, 'at hour ', hour, ' with output = ', output);
            }else if(output > 0 && edge.direction == 'negative'){
                fromnode = edge.to;
                tonode = edge.from;
                direction = "positive";
                console.log('Reverse Direction for Edge', id, 'at hour ', hour, ' with output = ', output);
            }
            edges.update({ id: id, width: normalizedWidth, from: fromnode, to: tonode, direction: direction, powerFlow: abs_output, label: Number(abs_output.toFixed(2)).toString()});
        }
    });


    network.redraw();
}

/*
document.getElementById('selectBox1').addEventListener('change', function() {
    console.log('Selected Option 1:', this.value);
});
document.getElementById('selectBox2').addEventListener('change', function() {
    console.log('Selected Option 2:', this.value);
});*/

// Call this when user load the page
function drawLegend() {

    var legendNodes = new vis.DataSet([
        {id: 1000, label: "Generator", group: "generator", value: 1},
        {id: 1001, label: "Wind Turbine", group: "windTurbine", value: 1},
        {id: 1002, label: "Solar Panel", group: "solarPanel", value: 1},
        {id: 1003, label: "Battery", group: "batteryStorage", value: 1},
        {id: 1004, label: "Load", group: "criticalLoad", value: 1},
        {id: 1005, label: "Node", group: "Node", value: 1}

    ]);
    var legendEdges = new vis.DataSet([]); // No edges for the legend
    var legendData = {
        nodes: legendNodes,
        edges: legendEdges
    };
    var legendOptions = {
        interaction: { dragNodes: false, dragView: false, zoomView: false }, // Disable all interactions
        layout: {
        hierarchical: {
            direction: 'LR',
            sortMethod: 'directed',
            treeSpacing: 70 // Adjust distance between trees
            }
        },
        nodes: {
            scaling: {
                min: 1,
                max: 1, // Set min and max scaling to 1 to prevent icon scaling
            },
        },
        physics: {
            enabled: false
        },
        groups: {
          generator: {
            shape: 'icon',
            icon: {
                face: 'Material Icons',
                code: '\ue932',
                size: 40,
                color: "#2B7CE9",
            }
          },
          windTurbine: {
            shape: 'icon',
            icon: {
                face: 'Material Icons',
                code: '\uec0c',
                size: 40,
                color: LIGHTBLUE,
            }
          },
          solarPanel: {
            shape: 'icon',
            icon: {
                face: 'Material Icons',
                code: '\uec0f',
                size: 40,
                color: YELLOW,
            }
          },
          batteryStorage: {
            shape: 'icon',
            icon: {
                face: 'Material Icons',
                code: '\ue1a3',
                size: 40,
                color: "#1cc88a",
            }
          },
          criticalLoad: {
            shape: 'icon',
            icon: {
                face: 'Material Icons',
                code: '\uea40',
                size: 40,
                color: "#4e73df",
            }
          },
          Node: {
            shape: 'icon',
            icon: {
                face: 'Material Icons',
                code: '\uef4a',
                size: 40,
                color: "#858796",
            }
          },
        },
    };

    var legendNetwork = new vis.Network(document.getElementById('legendNetwork'), legendData, legendOptions);

}