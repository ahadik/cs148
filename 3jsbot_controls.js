//////////////////////////////////////////////////
/////     ROBOT DYNAMICS PLACEHOLDER 
//////////////////////////////////////////////////

// CS148: add kinematics update from controls here
function robot_apply_controls() {
// apply robot controls to robot kinematics transforms and joint angles, then zero controls
// (for now) includes update of camera position based on base movement 

	robot.origin.xyz[0] += robot.control.xyz[0];
	robot.origin.xyz[1] += robot.control.xyz[1];
	robot.origin.xyz[2] += robot.control.xyz[2];
	robot.origin.rpy[0] += robot.control.rpy[0];
	robot.origin.rpy[1] += robot.control.rpy[1];
	robot.origin.rpy[2] += robot.control.rpy[2];
	
	robot.control.xyz = [0,0,0];
	robot.control.rpy = [0,0,0];
	
	
//	console.log(robot.joints[active_joint].servo.desired);
	//robot.joints[active_joint].servo.desired = robot.joints[active_joint].control;
	
	for (i in robot.joints){
		robot.joints[i].angle = robot.joints[i].angle + robot.joints[i].control;
		robot.joints[i].control = 0;
	}
	

    // move camera with robot base 
    // CS148: do not delete this
    // (need to pull this out and into 3jsbot support, at some point)
    camera_controls.object.position.x += robot.control.xyz[0];
    camera_controls.object.position.y += robot.control.xyz[1];
    camera_controls.object.position.z += robot.control.xyz[2];

}


