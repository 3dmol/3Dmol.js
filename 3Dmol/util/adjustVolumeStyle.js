import Gradient from "../Gradient";
import VolumeData from "../VolumeData";

// standardize voldata/volscheme in style
export default function adjustVolumeStyle(style){
    if (style) {
      if (style.volformat && !(style.voldata instanceof VolumeData)) {
        style.voldata = new VolumeData(style.voldata, style.volformat);
      }
      if (style.volscheme) {
        style.volscheme = Gradient.getGradient(style.volscheme);
      }
    }
  };