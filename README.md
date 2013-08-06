"Splinify" plugin for CINEMA 4D
=============

#[Download](https://github.com/eighteight/Splinify/archive/master.zip)
![Screenshot](https://raw.github.com/eighteight/Splinify/master/screenshot.png)
###[Video](https://vimeo.com) to follow

##Compability
The plugin has been compiled with the R14 sdk and has only been tested with C4D R14.

##Installation
Unzip the folder to your CINEMA/plugins/ folder and restart CINEMA.


##Usage
The plugin contains two objects that can be found under the ../plugins->Splinify menu entry.
###Splinify
A generator (place objects under it) that takes two objects and connects their points trying to minimize the overall distance.


##Settings
###Splinify
* Maximum Segment Length
  * Delete segments longer than this. Useful for getting rid of occasional lines connecting different parts of the model.
* Relative
  * When checked the _Maximum Segment Length_ is a multiplier of neighboring points, otherwise it's an absolute length.

