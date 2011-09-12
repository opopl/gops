
let angle_count=0

for ((x_angle=$xamin; x_angle<=$xamax; x_angle=$(($x_angle+$xai)) )); do	
for ((y_angle=$yamin; y_angle<=$yamax; y_angle=$(($y_angle+$yai)) )); do	
for ((z_angle=$zamin; z_angle<=$zamax; z_angle=$(($z_angle+$zai)) )); do	
					use exjmol
					angle_count=$(($angle_count+1))
done	
done		
done		

