# UI Tests 

UI testing is divided into two part. First **expect** and **action**, `expect` mean to look for the described output, if the desired expectations are not met then the ui testing is failed for the block in which `expect` is called. `action` refers to the interactions that user needs to do on the 3Dmol ui to reach the desired expectations. Inability of performing the described action means the ui testing has failed for the block in which action is called.

## Test 1 : Loading data from data base 

**expect** : The following ui must be present at the top of the viewer on call of `viewer.ui.initiateUI()` where viewer is created by using `3Dmol.createViewer` function call. 

![image-20210829160953969](Manual_test.assets/image-20210829160953969.png)

### Test 1.1 : Incorrect input submission 

**action** : click on submit button (green circular button with tick icon)

**expect** : Follwing change must be observed 

![image-20210829161127081](Manual_test.assets/image-20210829161127081.png)

## Test 2 : Loading files using model bar  

### Test 2.1 : New Model addition in 3Dmol js

**action** Change select to `pdb `and enter text `6zzm`in the input text area and press `enter` or submit button.

![image-20210829161832445](Manual_test.assets/image-20210829161832445.png)

**expect** Follwoing change in the title bar must be observed and UI must load a line figure of the whole model

![image-20210829161909400](Manual_test.assets/image-20210829161909400.png)

**action** : repeat the action mentioned above but change input fields as per the following table 

**expect** images in the third column must be expected 

| Drop down list value | Text input value | Desired Output                                               |
| -------------------- | ---------------- | ------------------------------------------------------------ |
| mmtf                 | 4hhb             | <img src="Manual_test.assets/image-20210829162356415.png" alt="image-20210829162356415" style="zoom:33%;" /> |
| cid                  | 2244             | <img src="Manual_test.assets/image-20210829162440501.png" alt="image-20210829162440501" style="zoom:33%;" /> |

 

### Test 2.2 : Hide/Show Model Toolbar

**action** : Click on the button with molecule as a symbol 

![image-20210829165547990](Manual_test.assets/image-20210829165547990.png)

**expect** : Title bar hidden 

![image-20210829165619049](Manual_test.assets/image-20210829165619049.png)

**action** : Click on the button with molecule as a symbol 

**expect** : Title bar shown

![image-20210829170205279](Manual_test.assets/image-20210829170205279.png)

## Test 3 : Selection Box

**action** : Load molecule 6zzm using model toolbar as shown in **Test 1**

**expect** : Molecule loaded and line structure is shown in the viewport 

![image-20210829200106265](Manual_test.assets/image-20210829200106265.png)

**action** :  Click on button with hand

**expect** : Plus button below it will be shown 

![image-20210829200547771](Manual_test.assets/image-20210829200547771.png)

### Test 3.1 : Addition of new selection before completion of previous form

**action** : Click on plus button 

**expect** : New selection card will be added 

![image-20210829200614059](Manual_test.assets/image-20210829200614059.png)

**action** : Click on plus button 

**expect** : Warning to complete the form should be shown

![image-20210829200651985](Manual_test.assets/image-20210829200651985.png)

### Test 3.2 : Addition of new selection : all atom 

**action** : Click on plus button 

**expect** : New selection card will be added 

**action** : Check Select All atom 

**expect** : All other option becomes invisible 

![image-20210829200732849](Manual_test.assets/image-20210829200732849.png)

**action** : Click on submit button 

**expect** : A new selection similar to follwoing image will be created with different id 

![image-20210829200752361](Manual_test.assets/image-20210829200752361.png)

### Test 3.3 : Addition of new selection : other selection 

**action** : Click on plus btton 

**expect** : New selection card will be created 

**action** : Check elem option 

**expect** : Input field wll be visible 

![image-20210829200859928](Manual_test.assets/image-20210829200859928.png)

#### Test 3.3.1 : Empty Selection 

**action** : Enter `NotAtom` as input 

**action** : click submit button 

**expect** : Error No atom selected 

![image-20210829200937608](Manual_test.assets/image-20210829200937608.png)

#### Test 3.3.2 : Carbon Atom only 

**action** : Enter `C` as input and press submit 

**expect** : New selection created with id 

![image-20210829201025479](Manual_test.assets/image-20210829201025479.png)

### Test 3.4 : Addition of new selection : Range Input 

**action** : Click on plus button 

**expect** : New selection card will be created 

**action** : Check `resi` and input range separated by comma `431-500,550-600`

![image-20210829201458270](Manual_test.assets/image-20210829201458270.png)

**action** : Click on submit 

**expect** : New Card will be created 

![image-20210829201549069](Manual_test.assets/image-20210829201549069.png)

### Test 3.5 : Addition of new selection : hetflag

**action** : Click on plus button 

**expect** : New selection card will be created 

**action** : Check `hetflag` and click submit 

![image-20210829201638637](Manual_test.assets/image-20210829201638637.png)

**expect** : New Card will be created 

![image-20210829201655148](Manual_test.assets/image-20210829201655148.png)

### Test 3.6 : Addition of new selection : Only carbon hetero atom 

**action** : Click on plus button 

**expect** : New selection card will be created 

**action** : Check `hetflag`, Check `elem` and enter `C` in the input field and click submit 

![image-20210829201748894](Manual_test.assets/image-20210829201748894.png)

**expect** : New Card will be created 

![image-20210829201815165](Manual_test.assets/image-20210829201815165.png)

### Test 3.7 : Addition of new selection : Chain A

**action** : Click on plus button 

**expect** : New selection card will be created 

**action** : Check `chain` and enter `A` in the input field and click submit 

![image-20210829201845380](Manual_test.assets/image-20210829201845380.png)

**expect** : New Card will be created 

![image-20210829201907435](Manual_test.assets/image-20210829201907435.png)

### Test 3.8 : Addition of new selection : Chain B

**action** : Click on plus button 

**expect** : New selection card will be created 

**action** : Check `chain` and enter `B` in the input field and click submit 

![image-20210829201944203](Manual_test.assets/image-20210829201944203.png)

**expect** : New Card will be created 

![image-20210829202008011](Manual_test.assets/image-20210829202008011.png)

## Test 4 : Addition of style 

### Test 4.1 : Addition of line style 

**action** : Click on plus atom in the selection card generated from Test 3.2 

![image-20210829202053602](Manual_test.assets/image-20210829202053602.png)

**action** : Check `line `and press submit 

![image-20210829202119550](Manual_test.assets/image-20210829202119550.png)

**expect** : Viewport has all the atoms drawn as line structure 

![image-20210829202152743](Manual_test.assets/image-20210829202152743.png)

### Test 4.2 : Addition of cross style 

**action** : Click on plus atom in the selection card generated from Test 3.5 

**action** : Check `cross `

**action** : Check `colorscheme` and select `ssPyMol `from the dropdown list 

![image-20210829202615337](Manual_test.assets/image-20210829202615337.png)

**expect** : Viewport has hetero atom drawn as cross structure 

![image-20210829202727455](Manual_test.assets/image-20210829202727455.png)

### Test 4.3 : Addition of stick style 

**action** : Click on plus atom in the selection card generated from Test 3.7 

**action** : Check `stick `

**action** : Check `colorscheme` and select `orangeCarbon`from the dropdown list 

**action** : Check `singleBonds` 

![image-20210829203232182](Manual_test.assets/image-20210829203232182.png)

**expect** : Viewport has chain A atom drawn as stick structure 

![image-20210829203307314](Manual_test.assets/image-20210829203307314.png)

### Test 4.4 : Addition of sphere style 

**action** : Click on plus atom in the selection card generated from Test 3.6 

**action** : Check `sphere `

**action** : Check `color` and select any color with red shade from the color menu 

**action** : Check `radius` and enter `0.4` and press submit

![image-20210829203358661](Manual_test.assets/image-20210829203358661.png)

**expect** : Viewport has hetero carbon atom drawn as red sphere structure 

![image-20210829203458076](Manual_test.assets/image-20210829203458076.png)

### Test 4.5 : Addition of Cartoon Style 

#### Test 4.5.1 : Cartoon Spectrum 

**action** : Click on plus atom in the selection card generated from Test 3.7 

**action** : Check `cartoon `

**action** : Check `color` and check `spectrum`

![image-20210829203605420](Manual_test.assets/image-20210829203605420.png)

**expect** : Viewport has chain A atom drawn as cartoon with spectrum colors. 

![image-20210829203629873](Manual_test.assets/image-20210829203629873.png)

#### Test 4.5.2 : Cartoon Colored 

**action** : Click on plus atom in the selection card generated from Test 3.8 

**action** : Check `cartoon `

**action** : Check `color` and select any color of your choice 

![image-20210829203850316](Manual_test.assets/image-20210829203850316.png)

**expect** : Viewport has chain B atom drawn as cartoon with selected color

![image-20210829203935329](Manual_test.assets/image-20210829203935329.png)

## Test 5 : Addition of new surface 

**action** : Click on surface button 

![image-20210829204217905](Manual_test.assets/image-20210829204217905.png)

**expect** : plus button shown 

![image-20210829204232344](Manual_test.assets/image-20210829204232344.png)

### Test 5.1 : Surface with default options 

**action** : Click on plus button 

**expect** : Surface card created 

![image-20210829204250336](Manual_test.assets/image-20210829204250336.png)

**action** : Click on plus button 

**expect** : Warning shown to complete the first card 

![image-20210829204309832](Manual_test.assets/image-20210829204309832.png)

**action** : Click submit button 

**expect** : VDW type surface added for all the atom and surface card with id created

![image-20210829204345367](Manual_test.assets/image-20210829204345367.png)

action : Click on the `-` symbol on the surface to remove this surface 

![image-20210829204459583](Manual_test.assets/image-20210829204459583.png)

expect : Surface created will be removed 

![image-20210829204521495](Manual_test.assets/image-20210829204521495.png)

### Test 5.2 : Surface with opacity for all atom 

**action** : Click on plus button 

**expect** : Surface card created 

**action** : Select opacity to some other value than 1 

**action** : Click submit button 

![image-20210829204628807](Manual_test.assets/image-20210829204628807.png)

**expect** : VDW type surface added for all the atom and surface card with id created

![image-20210829204731972](Manual_test.assets/image-20210829204731972.png)

### Test 5.3 : Surface made from self atom 

**action** : click on plus button 

**expect** : Surface card created 

**action** : Select `opacity `to around 0.6

**action** : Check `color` and set color to some shade of light blue 

**action** : Select self in surface atoms 

**action** : Select id created by Test 3.5 in show atom 

![image-20210829205134517](Manual_test.assets/image-20210829205134517.png)

**action** : Click `submit` option 

**expect** : VDW surface will be created around hetero carbon atoms

![image-20210829205157744](Manual_test.assets/image-20210829205157744.png)



## Test 6 :  Labels 

### Test 6.1 : Add atom label error message 

**action** : Click on spheres created by test 4.4 and add atom labels by clicking submit button 

![image-20210829205313884](Manual_test.assets/image-20210829205313884.png)

**expect** : Error message in the context menu 

![image-20210829205352363](Manual_test.assets/image-20210829205352363.png)

**action** : Click outside to hide the context menu 

![image-20210829205410849](Manual_test.assets/image-20210829205410849.png)

### Test 6.2 : Add atom label 

**action** : Click on any atom which can be clicked on the screen

**expect** : Context menu opened 

**action** : select property `x, y, z` and click submit button 

![image-20210829205447164](Manual_test.assets/image-20210829205447164.png)

**expect** : Label added at the atom location with `x, y, z` shown in the atom label 

![image-20210829205508995](Manual_test.assets/image-20210829205508995.png)

![image-20210829205525480](Manual_test.assets/image-20210829205525480.png)

### Test 6.3 : Remove atom label 

**action** : Click on atom with the label 

**expect** : Context menu opened with `remove atom label` option 

![image-20210829205551070](Manual_test.assets/image-20210829205551070.png)

**action** : Click on `Remove Label` option 

**expect** : Context menu hidden and label removed from the atom 

![image-20210829205639439](Manual_test.assets/image-20210829205639439.png)

### Test 6.4 : Add Selection Label 

**action** : Right click on any place in the viewport 

**expect** : Context Menu shows up 

**action** : Click `Add Label `

**action** : Input `Label Text`: Hetero Atom 

**action** : Set `Label Style` : Purple 

**action** : Set selection generated by Test 3.6

![image-20210829205838443](Manual_test.assets/image-20210829205838443.png)

**action** : Click submit button 

**expect** : A label added at the centroid of the atoms selected by selection of Test 3.6

![image-20210829205922894](Manual_test.assets/image-20210829205922894.png)

## Test 7 : Edit Selections and Style 

**action** : click on the button with pencil in the style card 

**expect** : parameters for the style are shown 

![image-20210829210001754](Manual_test.assets/image-20210829210001754.png)

**action** : Change the style value of any of the previously created style and click submit 

![image-20210829210042393](Manual_test.assets/image-20210829210042393.png)

**expect** : The changes will be reflected in the view port 

![image-20210829210116661](Manual_test.assets/image-20210829210116661.png)

**action** : Change the selection value of any of the previously created selection and click submit in the similar way

![image-20210829210139731](Manual_test.assets/image-20210829210139731.png)

![image-20210829210244664](Manual_test.assets/image-20210829210244664.png)

**expect** : The changes will be reflected in the view port in the styles drawn using that selection (no change on the previously generated surfaces)

![image-20210829210314389](Manual_test.assets/image-20210829210314389.png)

## Test 8 : Hide and remove style and selection 

**action** : Click on the eye button show on any of the selection 

![image-20210829210401999](Manual_test.assets/image-20210829210401999.png)

**expect** : Objects drawn through that selection will be hidden 

![image-20210829210503811](Manual_test.assets/image-20210829210503811.png)

![image-20210829210430859](Manual_test.assets/image-20210829210430859.png)

**action** : Click on the eye button shown on the selection 

![image-20210829210544446](Manual_test.assets/image-20210829210544446.png)

**expect** : Selection shown again 

![image-20210829210634051](Manual_test.assets/image-20210829210634051.png)

**action** : Repeat the same for any style in the selection 

![image-20210829210835132](Manual_test.assets/image-20210829210835132.png)

**expect** : Similar behaviour as shown above expect the ui should now only toggle a particular style 

![image-20210829210854562](Manual_test.assets/image-20210829210854562.png)

![image-20210829210928727](Manual_test.assets/image-20210829210928727.png)

**action** : Click on `-` button in the style card 

![image-20210829211006868](Manual_test.assets/image-20210829211006868.png)

**expect** : Style card and drawing will be removed 

![image-20210829211037105](Manual_test.assets/image-20210829211037105.png)

**action** : Click `-` button on the selection card

![image-20210829211107100](Manual_test.assets/image-20210829211107100.png)

**expect** : Selection card will be removed 

![image-20210829211141993](Manual_test.assets/image-20210829211141993.png)

## Test 9 : Edit surfaces 

**action** : Click the pencil icon on the surface card 

![image-20210829211229972](Manual_test.assets/image-20210829211229972.png)

**expect** : Parameters of the surface are shown 

![image-20210829211243387](Manual_test.assets/image-20210829211243387.png)

**action** : Change the surface parameters to some other acceptable parameters and click submit 

![image-20210829211455411](Manual_test.assets/image-20210829211455411.png)

**except** : New surface generated determined by the changed specification 

![image-20210829211516005](Manual_test.assets/image-20210829211516005.png)

## Test 10 : Removal of surfaces, styles and selections

**action** : Click on `-` button on the cards on the display 

**expect** : On clicking `-` card the figures generated by that card will be removed 

![image-20210829211553484](Manual_test.assets/image-20210829211553484.png)

![image-20210829211612014](Manual_test.assets/image-20210829211612014.png)



## Notes

* When the window is resized its overflow is updated but if the size of the window goes below 400 px then overflow padding that is coded in the lib does not work. For mobile devices the UI needs a redesign but works fine on the bigger screen

* The tests for data uri loading are saved in 3dmol-ui.html but are commented out, to use those tests uncomment the divs and compare the images for those divs 