//Event Handling
export class EventDispatcher {
  listeners = {} as Record<string, Array<(event: any) => void>>;

  dispatchEvent(event: any) {
    var listenerArray = this.listeners[event.type];

    if (listenerArray !== undefined) {
      event.target = this;

      // Disposing a resource removes its listener during this dispatch.
      const listeners = listenerArray.slice();
      for (var i = 0, l = listeners.length; i < l; i++)
        listeners[i].call(this, event);
    }
  }

  removeEventListener(type: string, listener?: (event: any) => void) {
    if (!listener) {
      this.listeners[type] = [];
    } else {
      const listeners = this.listeners[type];
      if (!listeners) return;
      var index = listeners.indexOf(listener);
      if (index !== -1) listeners.splice(index, 1);
    }
  }

  addEventListener(type: string, listener: (event: any) => void) {
    if (this.listeners[type] === undefined) this.listeners[type] = [];

    if (this.listeners[type].indexOf(listener) === -1)
      this.listeners[type].push(listener);
  }
}
