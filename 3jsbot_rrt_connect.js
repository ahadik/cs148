//////////////////////////////////////////////////
/////     RRT MOTION PLANNER
//////////////////////////////////////////////////

// CS148: 
// implement RRT-Connect by Kuffner and LaValle (2000)
//    paper link: http://msl.cs.uiuc.edu/~lavalle/papers/KufLav00.pdf

// compute motion plan and output into robot_path array 
// elements of robot_path are vertices based on tree structure in tree_init() 
// motion planner assumes collision checking by robot_collision_test()

/*
CS148: reference code has functions for:

    tree_add_vertex
    tree_add_edge
    random_config
    new_config
    nearest_neighbor
    rrt_extend
    rrt_connect
    find_path
    path_dfs
*/

var epsilon = .2;
var isConnected = false;
robot_path_traverse_idx = 0;

function colorPath(connectA, connectB){
  var workingNode = connectA;
  robot_path=[];
  
  while(workingNode.parent!=null){
    robot_path.unshift(workingNode);
    workingNode=workingNode.parent;
   // //console.log(workingNode);
  }

  workingNode = connectB;
  while(workingNode.parent!=null){
    robot_path.push(workingNode);
    workingNode=workingNode.parent;
      //  //console.log(workingNode);
  }

  draw_highlighted_path();
}

function draw_highlighted_path(){
		//console.log(robot_path);	
	for(var i=0;i<robot_path.length;i++){

		robot_path[i].geom.material.color = {r:1,g:0,b:0};
	}
	
}

//neighbor is a vertex object, desired is just the configuration points
function rrt_extend(tree, neighbor, desiredIn, connectBool){

  var desired;
  if(connectBool){
    desired = desiredIn.vertex;
  }else{
    desired = desiredIn;
  }

	var distance = getDistance(neighbor.vertex, desired);



	if((distance<epsilon)&&connectBool){

		isConnected = true;
    rrt_iterate=false;
		connect_trees(neighbor,desired);
    //console.log("CONNECTED");
    colorPath(neighbor,desiredIn);
		return false;
	}


	var xDiff = desired[0]-neighbor.vertex[0];
	var yDiff = desired[2]-neighbor.vertex[2];



	//Create a new vertex a distance of epsilon away from the neighbor

	var theta = Math.atan(yDiff/xDiff);

  var newY = Math.sin(theta)*epsilon+neighbor.vertex[2];
  var newX = Math.cos(theta)*epsilon+neighbor.vertex[0];


  if((yDiff<0)&&(xDiff<0)){
    newY = -Math.sin(theta)*epsilon+neighbor.vertex[2];
    newX = -Math.cos(theta)*epsilon+neighbor.vertex[0];
  }


  if((xDiff<0)&&(yDiff>0)){
    newY = -Math.sin(theta)*epsilon+neighbor.vertex[2];
    newX = -Math.cos(theta)*epsilon+neighbor.vertex[0];
  }
  

	var newVertex = [newX,0,newY];
	
	

	for(var i=3; i<neighbor.vertex.length;i++){
		
		newVertex.push(desired[i]);
	}

	//If the new vertex is a collision, return false
	if(robot_collision_test(newVertex)){
		    	////console.log(new Error().lineNumber);
		return false;
	}else{
		    	////console.log(new Error().lineNumber);
		tree_add_vertex(neighbor,newVertex, tree);
		return true;
	}
}

function getDistance(pointA, pointB){
	var dist = Math.sqrt(Math.pow((pointA[0]-pointB[0]),2)+Math.pow((pointA[2]-pointB[2]),2));

	return dist;
}

function nearest_neighbor(q,tree){
	////console.log(q);
	////console.log(tree);
	var nearest = [Infinity,Infinity];
	var nearestDist = Infinity;

	for (vertex in tree.vertices){
		var testDist = getDistance(q, tree.vertices[vertex].vertex);
		if(testDist<nearestDist){
			nearestDist = testDist;
			nearest = tree.vertices[vertex];
		}
	}
	return nearest;
}

function findRandom(){
	var xDist = Math.abs(robot_boundary[0][0])+Math.abs(robot_boundary[1][0]);
	var yDist = Math.abs(robot_boundary[0][2])+Math.abs(robot_boundary[1][2]);
	var x = (Math.random()*xDist)-Math.abs(robot_boundary[0][0]);
	var y = (Math.random()*yDist)-Math.abs(robot_boundary[0][2]);
	var r = 0;
	var p = Math.random()*7;
	var y = 0;
	

	var vertex = [x,0,y,r,p,y];
	for (x in robot.joints) {
        vertex.push(Math.random()*7);
    }

	
	return vertex;
}

function rrt_planning_iteration() {
  if(!isConnected){
    //////console.log("Not connected");
	//get a random vertex
	

	var vertex = findRandom();
	
	var treeA;
	var treeB;


	if(currTree==0){
		treeA = globalTrees[0];
		treeB = globalTrees[1];
		currTree=1;
	}else if(currTree==1){
		treeA = globalTrees[1];
		treeB = globalTrees[0];
		currTree=0;
	}
	//find the nearest neighbor on tree A to the randomly found vertex
	var nearestNeighbor = nearest_neighbor(vertex,treeA);


	//extend treeA towards the random vertex one step
	//Nearest neighbor is an object, vertex is just points
  //////console.log("single extend");
	rrt_extend(treeA, nearestNeighbor, vertex,false);
	////console.log(new Error().lineNumber);
	////console.log(treeA.vertices[treeA.newest].vertex);
	add_config_origin_indicator_geom(treeA.vertices[treeA.newest]);
	////console.log(new Error().lineNumber);
	//Find the vertex on treeB that is closest to the newly created vertex on treeB
	var nearestB = nearest_neighbor(treeA.vertices[treeA.newest].vertex,treeB);
////console.log(new Error().lineNumber);
	//While extension can proceed without collision, call rrt_extend extending treeA. If extend returns true, update nearestNeighbor to be the recently added element.

  var extendBool = true;

  var count = 0;
  while(extendBool){
    //////console.log("extending");

    if (rrt_iterate && (Date.now()-cur_time > 10)) {

        if(count == 10){break;}
        count++;
        extendBool = rrt_extend(treeB, nearestB, treeA.vertices[treeA.newest],true);
        if(isConnected){break;}
        if(extendBool){
          nearestB=treeB.vertices[treeB.newest];
          add_config_origin_indicator_geom(nearestB);
        }


        // update time marker for last iteration update
        cur_time = Date.now();
   }
   if(isConnected){break;}
	}
	if(isConnected){
		//connection has been made
    rrt_iterate=false;
		return true;
	}else{
		// connection has not been made
		return false;
	}
}
}

//parent is an object, child is just a vertex
function tree_add_vertex(parent,child, tree){
	////console.log(parent);
	////console.log(tree);
	//Add the new child to the tree
	tree.vertices[tree.newest+1] = {};
	tree.vertices[tree.newest+1].vertex = child;
	tree.vertices[tree.newest+1].edges = [];
	tree.vertices[tree.newest+1].parent = parent;
////console.log(new Error().lineNumber);
	//add the new child as a child of the parent
	parent.edges[parent.edges.length] = tree.vertices[tree.newest+1];
	tree.newest++;
}

function connect_trees(vertexA,vertexB){
	vertexA.connection = vertexB;
	vertexB.connection = vertexA;
}

function robot_rrt_planner_init() {

    // form configuration from base location and joint angles
    q_start_config = [
        robot.origin.xyz[0],
        robot.origin.xyz[1],
        robot.origin.xyz[2],
        robot.origin.rpy[0],
        robot.origin.rpy[1],
        robot.origin.rpy[2]
    ];

    q_names = {};  // store mapping between joint names and q DOFs

    for (x in robot.joints) {
        q_names[x] = q_start_config.length;
        q_start_config = q_start_config.concat(robot.joints[x].angle);
    }
    
    //console.log(q_start_config);

    // set goal configuration as the zero configuration
    var i; 
    q_goal_config = new Array(q_start_config.length);
    for (i=0;i<q_goal_config.length;i++) q_goal_config[i] = 0;

    // CS148: add necessary RRT initialization here
    globalTreeA = tree_init(q_start_config);
    globalTreeB = tree_init(q_goal_config);

    globalTrees = [globalTreeA, globalTreeB];
    currTree = 0;
    
    rrt_iterate = true;

    // make sure the rrt iterations are not running faster than animation update
    cur_time = Date.now();

    ////console.log("planner initialized");
}


function robot_rrt_planner_iterate() {

    rrt_alg = 1;  // 0: basic rrt (OPTIONAL), 1: rrt_connect (REQUIRED)

    if (rrt_iterate && (Date.now()-cur_time > 10)) {
        cur_time = Date.now();

        // CS148: implement RRT iteration here
        rrt_planning_iteration();

    }

    // return path not currently found
    return false;
}

function tree_init(q) {

    // create tree object
    var tree = {};

    // initialize with vertex for given configuration
    tree.vertices = [];
    tree.vertices[0] = {};
    tree.vertices[0].vertex = q;
    tree.vertices[0].edges = [];

    // create rendering geometry for base location of vertex configuration
    add_config_origin_indicator_geom(tree.vertices[0]);

    // maintain index of newest vertex added to tree
    tree.newest = 0;

    return tree;
}


function add_config_origin_indicator_geom(vertex) {
////console.log(new Error().lineNumber);
////console.log(vertex);
    // create a threejs rendering geometry for the base location of a configuration
    // assumes base origin location for configuration is first 3 elements 
    // assumes vertex is from tree and includes vertex field with configuration

    temp_geom = new THREE.CubeGeometry(0.1,0.1,0.1);
    temp_material = new THREE.MeshLambertMaterial( { color: 0xffff00, transparent: true, opacity: 0.7 } );
    temp_mesh = new THREE.Mesh(temp_geom, temp_material);
    ////console.log(new Error().lineNumber);
    temp_mesh.position.x = vertex.vertex[0];
    temp_mesh.position.y = vertex.vertex[1];
    temp_mesh.position.z = vertex.vertex[2];
    ////console.log(new Error().lineNumber);
    scene.add(temp_mesh);
////console.log(new Error().lineNumber);
    vertex.geom = temp_mesh;
}





