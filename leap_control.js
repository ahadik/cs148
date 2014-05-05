
		function vectorToString(vector){
			
			var vectorString = "";
			
			for(var i = 0; i<vector.length; i++){
				vectorString = vectorString+vector[i]+" ";
			}
			
			return vectorString;
		}
		
		var nullCoord = [.5,.5,.5];
		var palmCoord = undefined;
		
		function getCoord(){
			var returnVector;
			//If palmCoord has been set from the Leap loop, return its value, otherwise return null values of 0s.
			if (palmCoord!==undefined){
				returnVector=palmCoord;
			}else{
				returnVector=nullCoord;
			}
			return returnVector;
		}
	
		var controllerOptions = {enableGestures : true};
		Leap.loop(controllerOptions, function(frame){
		
			var interactionBox = frame.interactionBox;
		
			if(frame.hands.length > 0){
				var hand = frame.hands[0];
				
				palmCoord = interactionBox.normalizePoint(hand.palmPosition, true);
			}else{
				palmCoord = nullCoord;
			}	
		});
