"use strict"

//--------------vis network-------------
var nodes = null;
var edges = null;
var network = null;

var WIDTH_SCALE = 2,
  GREEN = "green",
  RED = "#C5000B",
  ORANGE = "orange",
  GRAY = "gray",
  BLACK = "#2B1B17",
  YELLOW = "#f6c23e",
  LIGHTBLUE = "#36b9cc";

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
  var legendNodes = [
    {id: 1000, x: x, y: y, label: "Generator", group: "generator", value: 1, fixed: true, physics: false},
    {id: 1001, x: x, y: y + step, label: "Wind Turbine", group: "windTurbine", value: 1, fixed: true, physics: false},
    {id: 1002, x: x, y: y + 2 * step, label: "Solar Panel", group: "solarPanel", value: 1, fixed: true, physics: false},
    {id: 1003, x: x, y: y + 3 * step, label: "Battery", group: "batteryStorage", value: 1, fixed: true, physics: false},
    {id: 1004, x: x, y: y + 4 * step, label: "Load", group: "criticalLoad", value: 1, fixed: true, physics: false},
    {id: 1005, x: x, y: y + 5 * step, label: "Node", group: "Node", value: 1, fixed: true, physics: false}
  ];
  // Add legend nodes to the nodes DataSet
  legendNodes.forEach(node => nodes.add(node));

// -----------------create a network-----------------
  var container = document.getElementById("mynetwork");
  var data = {
    nodes: nodes,
    edges: edges,
  };
  var options = {
    autoResize: true,
    nodes: {
      scaling: {
        min: 16,
        max: 32,
      },
    },
    edges: {
      color: GRAY,
      smooth: false,
      arrows: {
          to: { enabled: true, scaleFactor: 1, type: 'arrow' },
      }
    },
    physics: {
      barnesHut: { gravitationalConstant: -30000 },
      stabilization: { iterations: 2500 },
    },
    groups: {
      generator: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\ue932',
            size: 50,
            color: "#2B7CE9",
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
            color: "#1cc88a",
        }
      },
      criticalLoad: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\uec1b',
            size: 50,
            color: "#4e73df",
        }
      },
      Node: {
        shape: 'icon',
        icon: {
            face: 'Material Icons',
            code: '\uea0b',
            size: 50,
            color: "#858796",
        }
      },
    },
  };
  network = new vis.Network(container, data, options);
  // Fit the network once the stabilization is done
  // Function to check and adjust the view
    function checkAndAdjustView() {
        var scale = network.getScale();
        var viewPosition = network.getViewPosition();

        // Define your boundaries (example values)
        var minX = -400;
        var maxX = 100;
        var minY = -400;
        var maxY = 100;

        // Adjust scale if needed, for example, to prevent zooming out too far
        // Note: This is a simplistic approach. You may need a more complex logic based on your requirements
        if (scale < 0.05) { // Prevent zooming out too much
            network.moveTo({
                scale: 0.05
            });
        }

        // Adjust position if out of bounds
        if (viewPosition.x < minX || viewPosition.x > maxX || viewPosition.y < minY || viewPosition.y > maxY) {
            network.moveTo({
                position: {x: Math.min(Math.max(viewPosition.x, minX), maxX), y: Math.min(Math.max(viewPosition.y, minY), maxY)}
            });
        }
    }

    // Listen to dragEnd and zoom events to adjust the view
    network.on("dragEnd", function(params) {
        checkAndAdjustView();
    });

    network.on("zoom", function(params) {
        checkAndAdjustView();
    });

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
                           <p>Node Type: ${nodeData.group}</p>`;

            if (nodeData.group =="generator"){
                content += `<p>Generation Power: ${nodeData.value}W</p>`;
            }
            else{
                content += `<p>Nominal Voltage: ${nodeData.value}V</p>`;
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
                           <p>Length: ${edgeData.length}km</p>`;

            if (edgeData.edgetype == "Transformer"){
                content += `<p>Primary Voltage: ${edgeData.primaryVoltage}W</p>`;
                content += `<p>Secondary Voltage: ${edgeData.secondaryVoltage}W</p>`;
                content += `<p>Power Rating: ${edgeData.powerRating}W</p>`;
            }
            else if (edgeData.edgetype == "Overhead lines" || edgeData.edgetype == "Underground lines"){
                content += `<p>Power Flow: ${edgeData.powerFlow}W</p>`;
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




