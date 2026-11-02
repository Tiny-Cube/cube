use crate::{count, ensure, Edge, Graph, Model, Reaction, ReactionNetwork, Result};
// Deliberately restrict syntax to portable decimal tokens (no NaN, inf, hex, underscores).
pub fn integer(s: &str) -> Result<u64> {
    ensure(!s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()), "expected unsigned decimal integer")?;
    count(s.parse::<u64>().map_err(|_| "integer overflow")?)
}
pub fn number(s: &str) -> Result<f64> {
    ensure(!s.is_empty() && s.bytes().all(|b| b.is_ascii_digit() || b"+-.eE".contains(&b)), "expected decimal number")?;
    s.parse::<f64>().map_err(|_| "invalid decimal number".into())
}
struct Tokens<'a> { words: std::str::SplitWhitespace<'a> }
impl Tokens<'_> {
    fn next(&mut self) -> Result<&str> { self.words.next().ok_or("truncated KMC1 input".into()) }
    fn integer(&mut self) -> Result<u64> { integer(self.next()?) }
    fn number(&mut self) -> Result<f64> { number(self.next()?) }
}
pub fn loads(text: &str) -> Result<Box<dyn Model>> {
    let clean = text.lines().map(|l| l.split('#').next().unwrap_or("")).collect::<Vec<_>>().join(" ");
    let mut t = Tokens { words: clean.split_whitespace() };
    ensure(t.next()? == "KMC1", "expected KMC1 header")?;
    let kind = t.next()?.to_owned();
    let n = t.integer()?;
    let m = t.integer()?;
    let model: Box<dyn Model> = match kind.as_str() {
        "graph" => {
            let initial = t.integer()?;
            let mut edges = Vec::new();
            for _ in 0..m {
                edges.push(Edge { source: t.integer()?, target: t.integer()?, rate: t.number()?,
                    displacement: [t.number()?, t.number()?, t.number()?] });
            }
            Box::new(Graph::new(n, edges, initial)?)
        }
        "reaction" => {
            let mut populations = Vec::new();
            for _ in 0..n { populations.push(t.integer()?); }
            let mut reactions = Vec::new();
            for _ in 0..m {
                let rate = t.number()?;
                let mut reactants = Vec::new(); let mut products = Vec::new();
                for _ in 0..n { reactants.push(t.integer()?); }
                for _ in 0..n { products.push(t.integer()?); }
                reactions.push(Reaction { rate, reactants, products });
            }
            Box::new(ReactionNetwork::new(populations, reactions)?)
        }
        _ => return Err("unknown model kind".into()),
    };
    ensure(t.words.next().is_none(), "unexpected trailing tokens")?;
    Ok(model)
}
