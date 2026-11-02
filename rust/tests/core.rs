use universal_kmc::*;
use universal_kmc::events::{Event, EventModel};
struct Fixed(std::vec::IntoIter<f64>);
impl RandomSource for Fixed { fn uniform(&mut self) -> f64 { self.0.next().unwrap() } }
#[test]
fn clock_horizon_and_absorption() {
    let mut graph = Graph::new(1, vec![Edge{source:0,target:0,rate:2.0,displacement:[0.0;3]}, Edge{source:0,target:0,rate:3.0,displacement:[0.0;3]}],0).unwrap();
    let mut engine = Engine::with_rng(Fixed(vec![(-1.0_f64).exp(),0.8,(-1.0_f64).exp()].into_iter()));
    let mut rows = vec![];
    engine.run(&mut graph,10,0.3,|r| rows.push(r.clone())).unwrap();
    assert_eq!(rows.len(),3); assert_eq!(rows[1].event,1); assert_eq!(rows[2].event,-1); assert_eq!(engine.steps,1);
    assert_eq!(engine.time,0.3);
    let mut empty = Graph::new(1,vec![],0).unwrap();
    let mut engine = Engine::new(0).unwrap();
    assert!(engine.step(&mut empty,7.0).unwrap().is_none());
    assert_eq!(engine.time,7.0); assert_eq!(engine.status,"absorbed");
}
#[test]
fn reaction_conservation_and_validation() {
    let mut model = ReactionNetwork::new(vec![10,0],vec![Reaction{rate:0.1,reactants:vec![2,0],products:vec![0,1]},Reaction{rate:0.3,reactants:vec![0,1],products:vec![2,0]}]).unwrap();
    assert!((model.propensities().unwrap()[0]-4.5).abs()<1e-12);
    Engine::new(99).unwrap().run(&mut model,2000,f64::INFINITY,|r| assert_eq!(r.state[0]+2*r.state[1],10)).unwrap();
    assert!(Graph::new(1,vec![Edge{source:0,target:0,rate:-1.0,displacement:[0.0;3]}],0).is_err());
    assert!(io::loads("KMC1 graph 1 0 0 extra").is_err());
    assert!(rates::marcus(0.01,0.0,0.0,300.0).is_err());
}
#[test]
fn arbitrary_state_dynamic_catalog_and_minimal_model() {
    struct State { tasks: Vec<String>, value: f64 }
    let mut model = EventModel::new(State{tasks:vec!["a".into(),"b".into(),"c".into()],value:1.25}, |s: &State| {
        (0..s.tasks.len()).map(|i| Event::new(|_: &State| 1.0, move |s: &mut State| { s.tasks.remove(i); s.value*=0.5; })).collect()
    });
    let mut engine = Engine::new(13).unwrap();
    engine.run(&mut model,10,f64::INFINITY,|_|{}).unwrap();
    assert_eq!(engine.steps,3); assert_eq!(model.data.value,0.15625);
    struct Bare(u32);
    impl Model for Bare {
        fn propensities(&self)->Result<Vec<f64>> { Ok(if self.0<2 {vec![1.0]} else {vec![]}) }
        fn fire(&mut self,_:usize)->Result<()> {self.0+=1;Ok(())}
    }
    let mut bare=Bare(0);
    Engine::new(1).unwrap().run(&mut bare,10,f64::INFINITY,|_|{}).unwrap();
    assert_eq!(bare.0,2);
}
#[test]
fn rate_balance_and_rng() {
    let mut rng=PortableRng::new(0).unwrap();
    assert_eq!(rng.uniform(),48271.0/2147483647.0);
    let ratio=rates::marcus(0.01,0.2,0.05,300.0).unwrap()/rates::marcus(0.01,0.2,-0.05,300.0).unwrap();
    assert!((ratio-(-0.05/(rates::KB*300.0)).exp()).abs()<1e-13);
}
