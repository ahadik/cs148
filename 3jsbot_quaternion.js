//////////////////////////////////////////////////
/////     QUATERNION TRANSFORM ROUTINES 
//////////////////////////////////////////////////

/*
CS148: reference code has functions for:

quaternion_from_axisangle
quaternion_normalize
quaternion_multiply
quaternion_to_rotation_matrix
*/


function quaternion_from_axisangle(angle, axis){
	var ux = axis[0];
	var uy = axis[1];
	var uz = axis[2];
	
	var quaternion = [Math.cos(angle/2), ux*Math.sin(angle/2), uy*Math.sin(angle/2), uz*Math.sin(angle/2)];
	return quaternion;	
}

function quaternion_normalize(quaternion){
	var magnitude = Math.sqrt(quaternion[0]*quaternion[0]+quaternion[1]*quaternion[1]+quaternion[2]*quaternion[2]+quaternion[3]*quaternion[3]);
	
	var normalized_quaternion = [quaternion[0]/magnitude, quaternion[1]/magnitude, quaternion[2]/magnitude, quaternion[3]/magnitude];
	
	return normalized_quaternion;
}

function quaternion_to_rotation_matrix(q){
	var rotation_matrix = [
	[1-2*(q[2]*q[2]+q[3]*q[3]), 2*(q[1]*q[2]-q[0]*q[3]), 2*(q[0]*q[2]+q[1]*q[3]), 0],
	[2*(q[1]*q[2]+q[0]*q[3]), 1-2*(q[1]*q[1]+q[3]*q[3]), 2*(q[2]*q[3]-q[0]*q[1]), 0],
	[2*(q[1]*q[3]-q[0]*q[2]), 2*(q[0]*q[1]+q[2]*q[3]), 1-2*(q[1]*q[1]+q[2]*q[2]), 0],
	[0,0,0,1]
	];
	
	return rotation_matrix;
	
}