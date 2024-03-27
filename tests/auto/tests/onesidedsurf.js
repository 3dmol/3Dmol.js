$3Dmol.download('pdb:3erk', viewer, {}, function() {
  viewer.addSurface('MS',{'colorscheme':'whiteCarbon','onesided':true},
                {'byres':true, 'hetflag':false, 'within':{'distance':5,'sel':{'resn':'SB4'}}},
                {'hetflag':false});
  viewer.setStyle({'resn':'SB4'},{'stick':{'colorscheme':'chartreuseCarbon'}});
  viewer.setStyle({'bonds':0},{'sphere':{'radius':0.35}});
  viewer.zoomTo({'resn':'SB4'});
  viewer.render();
});
