"use strict"

// ---------Responsive-navbar-active-animation-----------
function test(){
	var tabsNewAnim = $('#navbarSupportedContent');
	var selectorNewAnim = $('#navbarSupportedContent').find('li').length;
	var activeItemNewAnim = tabsNewAnim.find('.active');
	var activeWidthNewAnimHeight = activeItemNewAnim.innerHeight();
	var activeWidthNewAnimWidth = activeItemNewAnim.innerWidth();
	var itemPosNewAnimTop = activeItemNewAnim.position();
	var itemPosNewAnimLeft = activeItemNewAnim.position();
	$(".hori-selector").css({
		"top":itemPosNewAnimTop.top + "px",
		"left":itemPosNewAnimLeft.left + "px",
		"height": activeWidthNewAnimHeight + "px",
		"width": activeWidthNewAnimWidth + "px"
	});
	$("#navbarSupportedContent").on("click","li",function(e){
		$('#navbarSupportedContent ul li').removeClass("active");
		$(this).addClass('active');
		var activeWidthNewAnimHeight = $(this).innerHeight();
		var activeWidthNewAnimWidth = $(this).innerWidth();
		var itemPosNewAnimTop = $(this).position();
		var itemPosNewAnimLeft = $(this).position();
		$(".hori-selector").css({
			"top":itemPosNewAnimTop.top + "px",
			"left":itemPosNewAnimLeft.left + "px",
			"height": activeWidthNewAnimHeight + "px",
			"width": activeWidthNewAnimWidth + "px"
		});
	});
}
$(document).ready(function(){
	setTimeout(function(){ test(); });
});
$(window).on('resize', function(){
	setTimeout(function(){ test(); }, 500);
});




// --------------add active class-on another-page move----------
jQuery(document).ready(function($){
	// Get current path and find target link
	var path = window.location.pathname.split("/").pop();

	// Account for home page with empty path
	if ( path == '' ) {
		path = 'main.html';
	}

	var target = $('#navbarSupportedContent ul li a[href="'+path+'"]');
	// Add active class to target link
	target.parent().addClass('active');
});


//--------------vis network-------------
var nodes = null;
var edges = null;
var network = null;

var WIDTH_SCALE = 2,
  GREEN = "green",
  RED = "#C5000B",
  ORANGE = "orange",
  GRAY = "gray",
  BLACK = "#2B1B17";

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
    {id: 1003, x: x, y: y + 3 * step, label: "Battery Storage", group: "batteryStorage", value: 1, fixed: true, physics: false},
    {id: 1004, x: x, y: y + 4 * step, label: "Critical Load", group: "criticalLoad", value: 1, fixed: true, physics: false},
    {id: 1005, x: x, y: y + 5 * step, label: "Controller", group: "controller", value: 1, fixed: true, physics: false}
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
        shape: "triangle",
        color: "#2B7CE9", // blue
      },
      windTurbine: {
        shape: "dot",
        color: "#5A1E5C", // purple
      },
      solarPanel: {
        shape: "square",
        color: "#C5000B", // red
      },
      batteryStorage: {
        shape: "square",
        color: "#FF9900", // orange
      },
      criticalLoad: {
        shape: "dot",
        color: "#109618", // green
      },
      controller: {
        shape: "dot",
        color: "#666666", // grey
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
    if (scale < 0.5) { // Prevent zooming out too much
        network.moveTo({
            scale: 0.5
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
}

//window.addEventListener("load", () => {
//  draw(microgridData);
//});


//andle the form submission to prevent the default action and make an AJAX request
$(document).ready(function() {
    $('#id_uploadGMLForm').submit(function(event) {
        event.preventDefault(); // Prevent the default form submission

        // Optional: Perform client-side validation here if needed

        // If the form is valid, make the AJAX call
        $.ajax({
            url: '/api/microgrid-data/', // The URL to your Django view
            type: 'GET', // or 'POST', depending on your view
            dataType: 'json',
            success: function(data) {
                console.log(data);
                // Process the data received from the server
                draw(data); // Assuming draw is your function to handle the data and display the graph
            },
            error: function(error) {
                console.log(error);
                // Handle errors here
            }
        });
    });
});



