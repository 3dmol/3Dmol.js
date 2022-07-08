//Event Handling
export class EventDispatcher {
  listeners = {};

  dispatchEvent(event) {
    var listenerArray = this.listeners[event.type];

    if (listenerArray !== undefined) {
      event.target = this;

      for (var i = 0, l = listenerArray.length; i < l; i++)
        listenerArray[i].call(this, event);
    }
  };

  removeEventListener(type, listener) {
    var index = this.listeners[type].indexOf(listener);

    if (index !== -1) this.listeners[type].splice(index, 1);
  };

  addEventListener(type, listener) {
    if (this.listeners[type] === undefined) this.listeners[type] = [];

    if (this.listeners[type].indexOf(listener) === -1)
      this.listeners[type].push(listener);
  };
}
