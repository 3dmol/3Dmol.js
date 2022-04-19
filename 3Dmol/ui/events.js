import $ from "jquery"

const Events = (function(){
  
    this.on = function(eventName, callback){
      $(document).on(eventName, (event)=>{
        callback(event);
      });
    }

    this.emit = function(eventName, eventData){

    }
  
});

export default Events;