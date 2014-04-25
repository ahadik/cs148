
function robot_set_planning_scene() {
    // currently just sets rendering geometries
    // world defined by robot_boundary and robot_obstacles objects in separate js include

    // set rendering geometries of world boundary
    temp_material = new THREE.MeshLambertMaterial( { color: 0xaf8c73, transparent: true, opacity: 0.4 } );

    temp_geom = new THREE.CubeGeometry(robot_boundary[1][0]-robot_boundary[0][0],0.2,0.2);
    temp_mesh = new THREE.Mesh(temp_geom, temp_material);
    temp_mesh.position.x = (robot_boundary[1][0]+robot_boundary[0][0])/2;
    temp_mesh.position.y = 0;
    temp_mesh.position.z = robot_boundary[0][2];
    scene.add(temp_mesh);

    temp_geom = new THREE.CubeGeometry(robot_boundary[1][0]-robot_boundary[0][0],0.2,0.2);
    temp_mesh = new THREE.Mesh(temp_geom, temp_material);
    temp_mesh.position.x = (robot_boundary[1][0]+robot_boundary[0][0])/2;
    temp_mesh.position.y = 0;
    temp_mesh.position.z = robot_boundary[1][2];
    scene.add(temp_mesh);

    temp_geom = new THREE.CubeGeometry(0.2,0.2,robot_boundary[1][2]-robot_boundary[0][2]);
    temp_mesh = new THREE.Mesh(temp_geom, temp_material);
    temp_mesh.position.x = robot_boundary[0][0];
    temp_mesh.position.y = 0;
    temp_mesh.position.z = (robot_boundary[1][2]+robot_boundary[0][2])/2;
    scene.add(temp_mesh);

    temp_geom = new THREE.CubeGeometry(0.2,0.2,robot_boundary[1][2]-robot_boundary[0][2]);
    temp_mesh = new THREE.Mesh(temp_geom, temp_material);
    temp_mesh.position.x = robot_boundary[1][0];
    temp_mesh.position.y = 0;
    temp_mesh.position.z = (robot_boundary[1][2]+robot_boundary[0][2])/2;
    scene.add(temp_mesh);
 
    // set rendering geometries of world obstacles
    var i;
    for (i=0;i<robot_obstacles.length;i++) { 
        temp_geom = new THREE.SphereGeometry(robot_obstacles[i].radius);
        temp_material = new THREE.MeshLambertMaterial( { color: 0xaf8c73, transparent: true, opacity: 0.6 } );
        temp_mesh = new THREE.Mesh(temp_geom, temp_material);
        temp_mesh.position.x = robot_obstacles[i].location[0][0];
        temp_mesh.position.y = robot_obstacles[i].location[1][0];
        temp_mesh.position.z = robot_obstacles[i].location[2][0];
        scene.add(temp_mesh);
    }
}


function robot_iscollision() {
    // test whether geometry of current configuration of robot is in collision with planning world 

    // form configuration from base location and joint angles
    var q_robot_config = [
        robot.origin.xyz[0],
        robot.origin.xyz[1],
        robot.origin.xyz[2],
        robot.origin.rpy[0],
        robot.origin.rpy[1],
        robot.origin.rpy[2]
    ];

    q_names = {};  // store mapping between joint names and q DOFs

    for (x in robot.joints) {
        	//console.log(new Error().lineNumber);
        q_names[x] = q_robot_config.length;
        q_robot_config = q_robot_config.concat(robot.joints[x].angle);
    }

    // test for collision and change base color based on the result
    if (robot_collision_test(q_robot_config)) {
        	//console.log(new Error().lineNumber);
       robot.links[robot.base].geom.material.color = {r:1,g:0,b:0};
       robot.links[robot.base].geom.material.opacity = 1.0;
    }
    else {
    	    	//console.log(new Error().lineNumber);
       robot.links[robot.base].geom.material.color = {r:0,g:0,b:1};
       robot.links[robot.base].geom.material.opacity = 0.7;
    }
         	//console.log(new Error().lineNumber);

}


function robot_collision_test(q) {
	//console.log("Collision");
	//console.log(q);

    // perform collision test of robot geometry against planning world 

    // test base origin (not extents) against world boundary extents
    if ((q[0]<robot_boundary[0][0])||(q[0]>robot_boundary[1][0])||(q[2]<robot_boundary[0][2])||(q[2]>robot_boundary[1][2])){
     	//console.log(new Error().lineNumber);
        return true;
    }
		    	//console.log(new Error().lineNumber);
    // traverse robot kinematics to test each body for collision
    return robot_collision_forward_kinematics(q);
}


function robot_collision_forward_kinematics (q) { 

    // transform robot base into the global world coordinates
        	//console.log(new Error().lineNumber);
    var mstack = matrix_multiply(generate_translation_matrix(q[0],q[1],q[2]), make_rotation_matrix(q[3],q[4],q[5]));
    //console.log("test");
    //var mstack = matrix_multiply(generate_translation_matrix(q[0],q[1],q[2]),matrix_multiply(matrix_multiply(generate_rotation_matrix_Z(q[5]),generate_rotation_matrix_Y(q[4])),generate_rotation_matrix_X(q[3])));
	//console.log(mstack);
    // recurse kinematics, testing collisions at each link
        	//console.log(new Error().lineNumber);
    return traverse_collision_forward_kinematics_link(robot.links[robot.base],mstack,q);
}


function traverse_collision_forward_kinematics_link(link,mstack,q) {
	//console.log(link);
	//console.log(mstack);
	//console.log(q);

    // test collision by transforming obstacles in world to link space
   // mstack_inv = matrix_invert_affine(mstack);
   
       	//console.log(new Error().lineNumber);
    mstack_inv = numeric.inv(mstack);
	//console.log("here");

    var i;
    var j;
    	//console.log(new Error().lineNumber);
    // test each obstacle against link bbox geometry by transforming obstacle into link frame and testing against axis aligned bounding box
    for (j=0;j<robot_obstacles.length;j++) { 
    	//console.log(new Error().lineNumber);
        var obstacle_local = matrix_multiply(mstack_inv,robot_obstacles[j].location);
    	//console.log(new Error().lineNumber);
        // assume link is in collision as default
        var in_collision = true; 
		
        // if obstacle lies outside the link extents along any dimension, no collision is detected
        if (
            (obstacle_local[0][0]<link.bbox.min.x-robot_obstacles[j].radius)
            ||
            (obstacle_local[0][0]>link.bbox.max.x+robot_obstacles[j].radius)
        ){
        		    	//console.log(new Error().lineNumber);
                in_collision = false;
        }
        if (
            (obstacle_local[1][0]<link.bbox.min.y-robot_obstacles[j].radius)
            ||
            (obstacle_local[1][0]>link.bbox.max.y+robot_obstacles[j].radius)
        )
                in_collision = false;
        if (
            (obstacle_local[2][0]<link.bbox.min.z-robot_obstacles[j].radius) 
            ||
            (obstacle_local[2][0]>link.bbox.max.z+robot_obstacles[j].radius)
        )
                in_collision = false;
        
        

        // if obstacle lies within link extents along all dimensions, a collision is detected and return true
        if (in_collision){
        	//console.log("1");
            return true;
        }
    }
    
   

    // recurse child joints for collisions, returning true if child returns collision
    if (typeof link.children !== 'undefined') { // return if there are no children
        	//console.log(new Error().lineNumber);
        for (i=0;i<link.children.length;i++) {
        	    	//console.log(new Error().lineNumber);
            if (traverse_collision_forward_kinematics_joint(robot.joints[link.children[i]],mstack,q)){
	            //console.log("2");
                return true;
            }
            
        }
    }

    // return false, when no collision detected for this link and children 
    //console.log("3");
    	//console.log(new Error().lineNumber);
    return false;
    //console.log("yada");
}


function traverse_collision_forward_kinematics_joint(joint,mstack,q) {
    	//console.log(new Error().lineNumber);
    var mstack_top = matrix_multiply(mstack,makeIdentityMatrix(4));
    
       	//console.log(new Error().lineNumber);

    // compute matrix transform origin of joint in the local space of the parent link
    var local_xform = matrix_multiply(generate_translation_matrix(joint.origin.xyz[0],joint.origin.xyz[1],joint.origin.xyz[2]), make_rotation_matrix(joint.origin.rpy[0],joint.origin.rpy[1],joint.origin.rpy[2]));
    
    //var local_xform = matrix_multiply(generate_translation_matrix(joint.origin.xyz[0],joint.origin.xyz[1],joint.origin.xyz[2]),matrix_multiply(matrix_multiply(generate_rotation_matrix_Z(joint.origin.rpy[2]),generate_rotation_matrix_Y(joint.origin.rpy[1])),generate_rotation_matrix_X(joint.origin.rpy[0])));
    	//console.log(new Error().lineNumber);
    // push local transform to origin to the top of the matrix stack
    var mstack_origin_top = matrix_multiply(mstack,local_xform)
    	//console.log(new Error().lineNumber);
    // transform motor rotation by quaternion for axis-angle joint rotation
    var tempvec = [joint.axis[0],joint.axis[1],joint.axis[2]]; 
        	//console.log(new Error().lineNumber);
    var tempquat = quaternion_from_axisangle(q[q_names[joint.name]],tempvec);
    var joint_local_quat = quaternion_normalize(tempquat);
    	//console.log(new Error().lineNumber);
    // push joint angle transform to the top of the matrix stack
    var joint_local_xform = quaternion_to_rotation_matrix(joint_local_quat);
        	//console.log(new Error().lineNumber);
    var mstack_top = matrix_multiply(mstack_origin_top,joint_local_xform); 

	//console.log("mstackTop");
	//console.log(mstack_top);
	//console.log(q);

    // recursively traverse child link with the current_xform being top of matrix stack 
    return traverse_collision_forward_kinematics_link(robot.links[joint.child],mstack_top,q);

}



