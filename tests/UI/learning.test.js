/** @jest-environment node */
jest.mock('canvas', () => ({}));
const { JSDOM } = require('jsdom');
const { LearningSession } = require('../../build/3Dmol.ui.js');
let dom, clients;
function fixture(config = {}, permissions = { edit: true, trajectory: true }) {
    const container = dom.window.document.createElement('div'); dom.window.document.body.append(container);
    const handlers = new Map(), labels = new Set(), atoms = [{ x: 0, y: 0, z: 0, chain: 'A', resn: 'ARG', resi: 12, atom: 'CA' }];
    const socket = { connected: true, on: jest.fn((event, fn) => handlers.set(event, fn)), off: jest.fn((event, fn) => { if (handlers.get(event) === fn) handlers.delete(event); }), emit: jest.fn(), disconnect: jest.fn() };
    const viewer = { container, selectedAtoms: () => atoms, getInternalState: () => ({models: [], labels: []}),
        getFrame: () => 0, getView: () => [0,0,0,0,0,0,0,1], getSlab: () => ({near:-10,far:10}),
        setClickable: jest.fn((sel, enabled, callback) => atoms.forEach(a => { a.clickable = enabled; if (callback) a.callback = callback; })),
        setHoverable: jest.fn((sel, enabled) => atoms.forEach(a => { a.hoverable = enabled; })),
        render: jest.fn(), setView: jest.fn(), getConfig:()=>({}), setFrame: jest.fn(async()=>{}), stopAnimate: jest.fn(),
        setInternalState: jest.fn(), subscribeStateChanges: jest.fn(()=>jest.fn()),
        addLabel: jest.fn((text, stylespec) => { const label={text,stylespec}; labels.add(label); return label; }),
        removeLabels: jest.fn(values => values.forEach(label=>labels.delete(label))) };
    const element = dom.window.document.createElement('div');container.append(element);
    const state = { viewer, ui: {element,getTrajectorySettings:()=>({playing:false}),setTrajectorySettings:jest.fn(async()=>{})},
        getViewSettings:()=>({background:'#ffffff'}), getControls:()=>({...permissions}),
        setControls:jest.fn(next=>{permissions={...permissions,...next};}), setDistancePicking:jest.fn(),
        subscribe:jest.fn(()=>jest.fn()),refresh:jest.fn(),setViewSettings:jest.fn(),
        modelTitles:new Map(),selections:new Map(),surfaces:new Map(),shapes:new Map(),volumes:new Map(),distances:new Map() };
    const session = new LearningSession(state, {...config,socketFactory:()=>socket}); clients.push(session);
    const receive=(event,value)=>handlers.get(event)?.(value);
    return {session,state,viewer,socket,container,atoms,labels,receive,handlers};
}
async function settle(){for(let i=0;i<8;i++)await Promise.resolve();}
beforeEach(()=>{jest.useFakeTimers();dom=new JSDOM('<!doctype html><body></body>',{url:'https://example.test/lesson'});clients=[];});
afterEach(()=>{clients.forEach(c=>c.destroy());dom.window.close();jest.useRealTimers();});

test('hidden mode requires an explicit room identifier',()=>{
    expect(()=>fixture({mode:'hidden'})).toThrow('requires a session identifier');
    expect(()=>fixture({mode:'hidden',session:'  '})).toThrow('requires a session identifier');
});

test('hidden viewers watch only their configured room and keep instructor/student roles independent',async()=>{
    const teacher=fixture({mode:'hidden',session:'lesson'}, {edit:false,trajectory:false});
    const student=fixture({mode:'hidden',session:'lesson'}), other=fixture({mode:'hidden',session:'other'});
    await Promise.all([teacher.session.ready,student.session.ready,other.session.ready]);
    for(const client of [teacher,student,other])expect(client.session.element.querySelector('.molui-learning-launch').hidden).toBe(true);
    expect(teacher.socket.emit).toHaveBeenCalledWith('watch session event',{name:'lesson'});
    expect(other.socket.emit).toHaveBeenCalledWith('watch session event',{name:'other'});
    teacher.container.focus();teacher.container.dispatchEvent(new dom.window.KeyboardEvent('keydown',{key:'i',ctrlKey:true,altKey:true,bubbles:true}));
    await settle();
    expect(teacher.socket.emit).toHaveBeenCalledWith('create session event',expect.objectContaining({name:'lesson'}));
    teacher.receive('session available',{name:'lesson'});teacher.receive('create session response',123);
    expect(teacher.session.getRole()).toBe('instructor');expect(teacher.state.getControls()).toEqual({edit:true,trajectory:true});
    student.receive('session available',{name:'lesson'});other.receive('session available',{name:'lesson'});await settle();
    expect(student.socket.emit).toHaveBeenCalledWith('join session event',{name:'lesson'});
    expect(other.socket.emit.mock.calls.some(([event])=>event==='join session event')).toBe(false);
    student.receive('join session response','student-id');teacher.receive('join session response','unsolicited');
    expect(teacher.session.getRole()).toBe('instructor');expect(student.session.getRole()).toBe('student');
    expect(student.state.getControls()).toEqual({edit:false,trajectory:false});expect(student.state.ui.element.hidden).toBe(true);
    expect(student.session.element.parentElement).toBe(student.container);expect(student.session.element.hidden).toBe(false);
    student.session.leave();expect(student.state.getControls()).toEqual({edit:true,trajectory:true});expect(student.state.ui.element.hidden).toBe(false);
    student.socket.emit.mockClear();student.receive('session available',{name:'lesson'});await settle();expect(student.socket.emit).not.toHaveBeenCalled();
    teacher.session.leave();expect(teacher.state.getControls()).toEqual({edit:false,trajectory:false});
});

test('hidden activation ignores text inputs, repeats and wrong shortcuts',async()=>{
    const c=fixture({mode:'hidden',session:'lesson'});await c.session.ready;c.socket.emit.mockClear();
    const input=c.session.element.querySelector('input');
    input.dispatchEvent(new dom.window.KeyboardEvent('keydown',{key:'i',ctrlKey:true,altKey:true,bubbles:true}));
    c.container.dispatchEvent(new dom.window.KeyboardEvent('keydown',{key:'i',ctrlKey:true,altKey:true,repeat:true,bubbles:true}));
    c.container.dispatchEvent(new dom.window.KeyboardEvent('keydown',{key:'i',ctrlKey:true,bubbles:true}));
    await settle();expect(c.socket.emit).not.toHaveBeenCalled();
});

test('losing simultaneous instructor activation joins the winning classroom',async()=>{
    const c=fixture({mode:'hidden',session:'lesson'});await c.session.start();
    c.receive('session available',{name:'lesson'});c.receive('create session response',0);await settle();
    expect(c.socket.emit).toHaveBeenCalledWith('join session event',{name:'lesson'});
    c.receive('join session response','student');expect(c.session.getRole()).toBe('student');
});

test.each(['start','join'])('disconnect during pending %s permits discovery and a fresh join',async(method)=>{
    const c=fixture({mode:'hidden',session:'lesson'});await c.session[method]();
    c.receive('disconnect','transport close');c.receive('connect');c.receive('session available',{name:'lesson'});await settle();
    c.receive('join session response','new-student');expect(c.session.getRole()).toBe('student');
});

test('server-ended hidden sessions keep watching for the next instructor',async()=>{
    const c=fixture({mode:'hidden',session:'lesson'});await c.session.join();c.receive('join session response','first');
    c.receive('leave session response',{reason:'instructor disconnected'});expect(c.session.getRole()).toBe('idle');
    c.receive('session available',{name:'lesson'});await settle();c.receive('join session response','second');expect(c.session.getRole()).toBe('student');
});

test('button mode prefills a name without joining or creating and supports normal controls',async()=>{
    const c=fixture({session:'lesson'});await c.session.ready;
    expect(c.socket.emit).not.toHaveBeenCalled();expect(c.session.element.querySelector('input').value).toBe('lesson');
    c.session.element.querySelector('.molui-learning-launch').click();
    expect(c.socket.emit).toHaveBeenCalledWith('check session name event',{name:'lesson'});
    c.receive('check session name response','free');
    const create=[...c.session.element.querySelectorAll('button')].find(b=>b.textContent==='Create');create.click();await settle();
    c.receive('create session response',123);expect(c.session.getRole()).toBe('instructor');
    const ask=[...c.session.element.querySelectorAll('button')].find(b=>b.textContent==='Query Atoms');ask.click();
    expect(c.socket.emit).toHaveBeenCalledWith('query start',{name:'lesson',kind:'atoms'});
    c.receive('query start response','atoms');ask.click();expect(c.socket.emit).toHaveBeenCalledWith('query fetch',{name:'lesson'});
    c.receive('query fetch response',[[[0,0,0],2]]);expect(c.labels.size).toBe(1);
    ask.click();expect(c.socket.emit).toHaveBeenCalledWith('query end',{name:'lesson'});
    c.receive('query end response');expect(c.labels.size).toBe(0);
});

test('student query callbacks cannot create answers after the query ends or the session leaves',async()=>{
    const c=fixture();await c.session.join('lesson');c.receive('join session response','student');
    c.receive('query start response','atoms');const callback=c.atoms[0].callback;callback(c.atoms[0]);
    expect(c.socket.emit).toHaveBeenCalledWith('query update',{name:'lesson',selected:[[0,0,0]]});expect(c.labels.size).toBe(1);
    c.receive('query end response');expect(c.atoms[0].clickable).toBe(false);expect(c.labels.size).toBe(0);
    c.socket.emit.mockClear();callback(c.atoms[0]);expect(c.socket.emit).not.toHaveBeenCalled();
    c.receive('query start response','atoms');const next=c.atoms[0].callback;c.session.leave();c.socket.emit.mockClear();next(c.atoms[0]);expect(c.socket.emit).not.toHaveBeenCalled();
});

test('destroy removes its own UI, transport handlers, shortcuts and restores original atom callbacks',async()=>{
    const c=fixture({mode:'hidden',session:'lesson'}), callback=jest.fn();
    c.atoms[0].clickable=true;c.atoms[0].callback=callback;
    await c.session.join();c.receive('join session response','student');c.receive('query start response','atoms');
    c.session.destroy();expect(c.atoms[0].callback).toBe(callback);expect(c.atoms[0].clickable).toBe(true);
    expect(c.container.querySelector('.molui-learning')).toBeNull();expect(c.socket.disconnect).toHaveBeenCalledTimes(1);expect(c.handlers.size).toBe(0);
    c.socket.emit.mockClear();c.container.dispatchEvent(new dom.window.KeyboardEvent('keydown',{key:'i',ctrlKey:true,altKey:true,bubbles:true}));await settle();expect(c.socket.emit).not.toHaveBeenCalled();
    c.session.destroy();expect(c.socket.disconnect).toHaveBeenCalledTimes(1);
});


test('hidden connection failures stay unobtrusive until the instructor activates the viewer',async()=>{
    const c=fixture({mode:'hidden',session:'lesson'});await c.session.ready;
    c.receive('connect_error',new Error('Server unavailable'));
    const dialog=c.session.element.querySelector('.molui-learning-dialog');
    expect(dialog.hidden).toBe(true);
    expect(c.session.element.querySelector('.molui-learning-launch').hidden).toBe(true);
    c.container.dispatchEvent(new dom.window.KeyboardEvent('keydown',{key:'i',ctrlKey:true,altKey:true,bubbles:true}));
    await settle();expect(dialog.hidden).toBe(false);expect(dialog.textContent).toContain('Server unavailable');
    c.session.leave();c.receive('connect_error',new Error('Still unavailable'));expect(dialog.hidden).toBe(true);
    c.session.open();expect(dialog.hidden).toBe(false);expect(dialog.textContent).toContain('Still unavailable');
});
