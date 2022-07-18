//Event Handling
export class EventDispatcher {
  listeners = {} as Record<string, Array<(event: any) => void>>;

  dispatchEvent(event: any) {
    var listenerArray = this.listeners[event.type];

    if (listenerArray !== undefined) {
      event.target = this;

      for (var i = 0, l = listenerArray.length; i < l; i++)
        listenerArray[i].call(this, event);
    }
  }

  removeEventListener(type: string, listener?: (event: any) => void) {
    if (!listener) {
      this.listeners[type] = [];
    } else {
      var index = this.listeners[type].indexOf(listener);
      if (index !== -1) this.listeners[type].splice(index, 1);
    }
  }

  addEventListener(type: string, listener: (event: any) => void) {
    if (this.listeners[type] === undefined) this.listeners[type] = [];

    if (this.listeners[type].indexOf(listener) === -1)
      this.listeners[type].push(listener);
  }
}
