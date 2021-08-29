# UI Tests 

UI testing is divided into two part. First **expect** and **action**, `expect` mean to look for the described output, if the desired expectations are not met then the ui testing is failed for the block in which `expect` is called. `action` refers to the interactions that user needs to do on the 3Dmol ui to reach the desired expectations. Inability of performing the described action means the ui testing has failed for the block in which action is called.

## Test 1 : Loading data from data base 

**expect** : The following ui must be present at the top of the viewer on call of `viewer.ui.initiateUI()` where viewer is created by using `3Dmol.createViewer` function call. 

![image-20210829160953969](Manual_test.assets/image-20210829160953969.png)

### Incorrect input submission 

**action** : click on submit button (green circular button with tick icon)

**expect** : Follwing change must be observed 

![image-20210829161127081](Manual_test.assets/image-20210829161127081.png)

## Loading files using model bar  

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

 

### Hide/Show Model Toolbar

**action** : Click on the button with molecule as a symbol 

![image-20210829165547990](Manual_test.assets/image-20210829165547990.png)

**expect** : Title bar hidden 

![image-20210829165619049](Manual_test.assets/image-20210829165619049.png)

**action** : Click on the button with molecule as a symbol 

**expect** : Title bar shown

![image-20210829170205279](Manual_test.assets/image-20210829170205279.png)

## Test 2 : Selection Box

**action** : Load molecule 6zzm using model toolbar as shown in **Test 1**

**expect** : Molecule loaded and line structure is shown in the viewport 

**action** :  Click on button with hand

**expect** : Plus button below it will be shown 

### Addition of new selection before completion of previous form

action : Click on plus button 

expect : New selection card will be added 

action : Click on plus button 

expect : Warning to complete the form should be shown

### Addition of new selection : all atom 

action : Click on plus button 

expect : New selection card will be added 

action : Check Selection All atom 

expect : All other option becomes invisible 

action : Click on submit button 

expect : A new selection similar to follwoing image will be created with different id 

### Addition of new selection : other selection 

action : Click on plus btton 

expect : New selection card will be created 

action : Check elem option 

expect : Input field wll be visible 

#### Empty Selection 

action : Enter `NotAtom`as input 

action : click submit button 

expect : Error invalid selection 

#### Carbon Atom only 

action : Enter `C` as input 

expect : New selection created with id 

### Addition of new selection : Range Input 

action : Click on plus btton 

expect : New selection card will be created 

action : 