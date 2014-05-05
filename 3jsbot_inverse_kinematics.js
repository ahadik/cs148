//////////////////////////////////////////////////
/////     INVERSE KINEMATICS 
/////     Resolved-rate IK with geometric jacobian
//////////////////////////////////////////////////

// CS148: generate joint controls to move robot to move robot endeffector to target location

/*
CS148: reference code has functions for:

robot_inverse_kinematics
iterate_inverse_kinematics
*/


function make_thing(matrix){
	var betterMax =[
		[1, 0 , 0, matrix[0][0]],
		[0, 1 , 0, matrix[1][0]],
		[0, 0 , 1, matrix[2][0]],
		[0, 0 , 0, 1],
	]
	
	return betterMax;
}

function robot_inverse_kinematics(target_pos, endeffector_joint, endeffector_local_pos) {
    // compute joint angle controls to move location on specified link to Cartesian location
    if (update_ik) {
    
    	var endeffector_mat = make_thing(multiply_matrices(endeffector_joint.xform,endeffector_local_pos));
    	simpleApplyMatrix(endeffector_geom, matrix_2Darray_to_threejs(endeffector_mat));
        
        
        simpleApplyMatrix(target_geom, matrix_2Darray_to_threejs([[1,0,0,target_pos[0][0]],[0,1,0,target_pos[1][0]],[0,0,1,target_pos[2][0]],[0,0,0,target_pos[3][0]]]));
		
		
		
		endeffector_geom.visible = true;
        target_geom.visible = true;
       	
       	iterate_inverse_kinematics(target_pos, endeffector_joint, endeffector_local_pos);
        
        

    }
    else {
        endeffector_geom.visible = false;
        target_geom.visible = false;
    }
    update_ik = false;

}

function getOn(joint, endeffector_global_pos){

	
	var end_x = endeffector_global_pos[0];
	var end_y = endeffector_global_pos[1];
	var end_z = endeffector_global_pos[2];

	var o_sub = [end_x-joint[0][3],end_y-joint[1][3],end_z-joint[2][3]];

	
	return o_sub;
}

function getJacobianVector (zi, o_sub){

	var vector_cros = vector_cross(zi,o_sub);

	var jacobianVector = [[vector_cros[0]],[vector_cros[1]],[vector_cros[2]],[zi[0]],[zi[1]],[zi[2]]];
	return jacobianVector;
}

function pseudo_inverse(jacobian, target_pos, global_end){
	alpha=0.1;
	var jacobianTrans = column_matrix_transpose(jacobian);
	var jacobianRow = column_to_row(jacobian);
	
	var dx = new Array(6);
	dx[0]=[target_pos[0][0]-global_end[0][0]];
	dx[1]=[target_pos[1][0]-global_end[1][0]];
	dx[2]=[target_pos[2][0]-global_end[2][0]];
	dx[3]=[0];
	dx[4]=[0];
	dx[5]=[0];
	
	pseudoInverseMult = matrix_multiply(jacobianTrans, jacobianRow);
	
	var pseudoInverse = numeric.inv(pseudoInverseMult);
	
	
	var jacobiansMult = matrix_multiply(pseudoInverse, jacobianTrans);
	var matrix = numeric.dot(jacobiansMult, dx);
	var finalMat = scalar_mult(matrix,alpha);
	
	return finalMat;
	
}


function jacobian_transpose(jacobian, target_pos, global_end){

	var alpha = 0.1;
	var jacobianTrans = column_matrix_transpose(jacobian);
	
	dx = new Array(6);
	dx[0]=[target_pos[0][0]-global_end[0][0]];
	dx[1]=[target_pos[1][0]-global_end[1][0]];
	dx[2]=[target_pos[2][0]-global_end[2][0]];
	dx[3]=[0];
	dx[4]=[0];
	dx[5]=[0];
	
		
	var pseudo = multiply_matrices(jacobianTrans,dx);

	
	var pseudo_scaled = scalar_mult(pseudo,alpha);
	
	return pseudo_scaled;

}

function getZi(joint){
	var joint = robot.joints[joint];
	var joint_xform = joint.xform;
	
	var zeroed_xform = new Array(4);
	
	zeroed_xform[0] = new Array(4);
	zeroed_xform[1] = new Array(4);
	zeroed_xform[2] = new Array(4);
	zeroed_xform[3] = new Array(4);
	
	zeroed_xform[0][0] = joint_xform[0][0];
	zeroed_xform[0][1] = joint_xform[0][1];
	zeroed_xform[0][2] = joint_xform[0][2];
	zeroed_xform[0][3] = 0;
	
	zeroed_xform[1][0] = joint_xform[1][0];
	zeroed_xform[1][1] = joint_xform[1][1];
	zeroed_xform[1][2] = joint_xform[1][2];
	zeroed_xform[1][3] = 0;
	
	zeroed_xform[2][0] = joint_xform[2][0];
	zeroed_xform[2][1] = joint_xform[2][1];
	zeroed_xform[2][2] = joint_xform[2][2];
	zeroed_xform[2][3] = 0;
	
	zeroed_xform[3][0] = joint_xform[3][0];
	zeroed_xform[3][1] = joint_xform[3][1];
	zeroed_xform[3][2] = joint_xform[3][2];
	zeroed_xform[3][3] = joint_xform[3][3];
	
	
	
	var mult_axis = [[joint.axis[0]],[joint.axis[1]],[joint.axis[2]],[1]];

	var zi = multiply_matrices(zeroed_xform,mult_axis);

	
	var zi_return = [zi[0][0],zi[1][0],zi[2][0]];
	
	return zi_return;
}

function applyPseudo(pseudo, joint){
	var workingJoint = joint;

	
	var count = 0;
	
	while(workingJoint.parent!="base"){
		workingJoint.control+=pseudo[count][0];

		count++;
		if(workingJoint.parent!=robot.base){
			workingJoint=robot.joints[robot.links[workingJoint.parent].parent];
		}
	}
	workingJoint.control+=pseudo[count][0];

	
}

function iterate_inverse_kinematics(target_pos, endeffector_joint, endeffector_local_pos){

		endeffector_global_pos_2 = multiply_matrices(endeffector_joint.xform, [[1,0,0,endeffector_local_pos[0]],[0,1,0,endeffector_local_pos[1]],[0,0,1,endeffector_local_pos[2]],[0,0,0,endeffector_local_pos[3]]]);
		
		endeffector_global_pos = multiply_matrices(endeffector_joint.xform, endeffector_local_pos);
        


	var jacobian = new Array(4);
	var jacobianVec;
		
	var workingJoint = endeffector_joint;
	var count = 0;
	while (workingJoint.parent!="base"){
		var zi = getZi(workingJoint.name);
		
		var o_sub = getOn(workingJoint.xform, endeffector_global_pos);
				
		jacobianVec = getJacobianVector (zi, o_sub);
		jacobian[count]=jacobianVec;
		workingJoint = robot.joints[robot.links[workingJoint.parent].parent];
		count++;
	}
	
	zi = getZi(workingJoint.name);

	o_sub = getOn(workingJoint.xform, endeffector_global_pos);
	jacobianVec = getJacobianVector (zi, o_sub);
	count++;		
	jacobian[3]=jacobianVec;


	
	//var deltatheta = jacobian_transpose(jacobian, target_pos, endeffector_global_pos);
	var deltatheta = pseudo_inverse(jacobian, target_pos, endeffector_global_pos);

	applyPseudo(deltatheta,endeffector_joint);
	//robot_pd_control();

	
}