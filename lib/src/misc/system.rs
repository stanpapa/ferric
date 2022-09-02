use std::error;
use std::process::Command;

// todo: create enum for binaries

pub fn system_command(s: &str) -> Result<(), Box<dyn error::Error>> {
    let command = format!("{}", s);

    let output = Command::new(command)
        .output()
        .unwrap_or_else(|e| panic!("failed to execute process: {}", e));

    if output.status.success() {
        let out = String::from_utf8_lossy(&output.stdout);
        print!("{}", out);
    } else {
        let out = String::from_utf8_lossy(&output.stderr);
        print!("{} failed and stderr was:\n{}", s, out);
    }

    Ok(())
}
