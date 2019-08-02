#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = Camera_SavingProcs


Static function go()
	wave w1=cam_0

	make/o/n=(dimsize(w1,0),dimsize(w1,1),dimsize(w1,2),1) Test
	variable n=0
	string name
	do
		sprintf name "Image%g",n
		wave w1=$name
		print nameofwave(w1)
		Test[][][][n]=w1[p][q][r]
		n+=1
		killwaves w1
		Insertpoints/M=3 (n), 1, Test




	while(n<123)

end

Static function gosavecolor()

	variable n=0
	string name
	do
		sprintf name "Image%g",n
		wave w1=$name
		imagetransform compress w1
		print nameofwave(w1)

		imagetransform compress w1
		wave W_Compressed
		Save/C/P=Hold1 W_Compressed as nameofwave(w1)+".ibw"
		n+=1
		//killwaves w1

	while(n<124)

end

static function gosavegrey()

	variable n=0
	string name
	do
		sprintf name "Image%g",n
		wave w1=$name
		imagetransform compress w1
		print nameofwave(w1)
		imagetransform rgb2gray w1
		wave M_RGB2Gray
		imagetransform compress M_RGB2Gray
		wave W_Compressed
		Save/C/P=Hold1 W_Compressed as nameofwave(w1)+".ibw"
		n+=1
		//killwaves w1

	while(n<123)

end

static function gosaveJPG()

	variable n=0
	string name
	newpath/q/o SavePath "D:Devin:Documents:Data:Boulder:SuperSmall:20160829_Camera:Scan1_Raw:JPG"
	do
		sprintf name "Image%g",n
		wave w1=$name
		duplicate/o w1 Image
		•Display/N=SaveWindow;AppendImage Image
		•ModifyImage Image ctab= {*,*,Rainbow,0}
•ModifyGraph margin(left)=29,margin(bottom)=29,margin(top)=14,margin(right)=14;DelayUpdate
•ModifyGraph width=288,height=216
•ModifyGraph tick=2,noLabel=2
•ModifyGraph margin=14
•ModifyGraph mirror=1
•ModifyGraph tick=3,mirror=0
SavePICT/o/p=SavePath/E=-6/B=72 as Name+".jpg"
killwindow SaveWindow
		//Save/C/P=Hold1 W_Compressed as nameofwave(w1)+".ibw"
		n+=1
		//killwaves w1

	while(n<125)

end
